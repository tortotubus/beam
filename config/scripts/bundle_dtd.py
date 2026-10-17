#!/usr/bin/env python3
"""Bundle local SYSTEM parameter-entity includes; not a general DTD processor."""
import argparse
from pathlib import Path
import re
import sys
from urllib.parse import unquote, urlsplit

NAME = r"[A-Za-z_:][A-Za-z0-9_.:-]*"
TOKEN = re.compile(r"\s+|<!--(?:[^-]|-[^-])*-->|<![^>\"']*(?:(?:\"[^\"]*\"|'[^']*')[^>\"']*)*>|%" + NAME + r";", re.DOTALL)
INCLUDE = re.compile(r"<!ENTITY\s+%\s+(" + NAME + r")\s+SYSTEM\s+([\"'])(.*?)\2\s*>", re.DOTALL)
ELEMENT = re.compile(r"<!ELEMENT\s+(" + NAME + r")\s+")
ATTLIST = re.compile(r"<!ATTLIST\s+" + NAME + r"\s+")


class BundleError(ValueError):
    pass


class Bundler:
    def __init__(self):
        self.entities = {}
        self.elements = {}
        self.dependencies = []
        self.stack = []

    def expand(self, path):
        path = Path(path).resolve()
        if path in self.stack:
            raise BundleError("Include cycle: " + " -> ".join(map(str, self.stack + [path])))
        if len(self.stack) >= 64:
            raise BundleError(f"{path}: include nesting exceeds 64 files")
        try:
            source = path.read_text(encoding="utf-8")
        except (OSError, UnicodeError) as error:
            raise BundleError(f"Cannot read DTD {path}: {error}") from error
        if path not in self.dependencies:
            self.dependencies.append(path)
        self.stack.append(path)
        # Preserve provenance without putting machine-specific absolute paths in artifacts.
        label = path.name.replace("--", "- -")
        output = [f"<!-- Begin {label} -->\n"]
        position = 0
        while position < len(source):
            location = f"{path}:{source.count(chr(10), 0, position) + 1}"
            match = TOKEN.match(source, position)
            if not match:
                raise BundleError(f"{location}: malformed or unsupported DTD syntax")
            token = match.group()
            position = match.end()
            if token.isspace() or token.startswith("<!--"):
                output.append(token)
                continue
            include = INCLUDE.fullmatch(token)
            if include:
                name, _, system = include.groups()
                if name in self.entities:
                    raise BundleError(f"{location}: duplicate parameter entity {name}; first at {self.entities[name][1]}")
                uri = urlsplit(system)
                if uri.scheme or uri.netloc or uri.query or uri.fragment or not uri.path:
                    raise BundleError(f"{location}: only local SYSTEM paths are supported: {system}")
                target = Path(unquote(uri.path))
                if target.is_absolute():
                    raise BundleError(f"{location}: SYSTEM paths must be relative: {system}")
                self.entities[name] = (path.parent / target, location)
            elif token.startswith("%"):
                name = token[1:-1]
                if name not in self.entities:
                    raise BundleError(f"{location}: undefined parameter entity {name}")
                target, declaration = self.entities[name]
                try:
                    output.append(self.expand(target))
                except BundleError as error:
                    raise BundleError(f"{location}: expanding %{name}; declared at {declaration}\n{error}") from error
            else:
                element = ELEMENT.match(token)
                if not element and not ATTLIST.match(token):
                    raise BundleError(f"{location}: unsupported declaration; use ELEMENT, ATTLIST or local SYSTEM parameter entities")
                if re.search(r"%" + NAME + r";", token):
                    raise BundleError(f"{location}: parameter entities must include complete declarations")
                if element:
                    name = element[1]
                    if name in self.elements:
                        raise BundleError(f"{location}: duplicate element {name}; first at {self.elements[name]}")
                    self.elements[name] = location
                output.append(token)
        self.stack.pop()
        output.append(f"\n<!-- End {label} -->\n")
        return "".join(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    try:
        bundler = Bundler()
        result = bundler.expand(args.root)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != result:
            args.output.write_text(result, encoding="utf-8")
        # CMake uses these dependencies to reconfigure when any include changes.
        for dependency in bundler.dependencies:
            print(dependency)
    except (BundleError, OSError) as error:
        print(error, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
