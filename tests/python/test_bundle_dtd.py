import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile
import time
import unittest

sys.dont_write_bytecode = True

SOURCE = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("bundle_dtd", SOURCE / "config/scripts/bundle_dtd.py")
bundle = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bundle)


class BundleDTDTest(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.root = Path(self.temporary.name)

    def write(self, name, content):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        return path

    def expand(self, text):
        return bundle.Bundler().expand(self.write("root.dtd", text))

    def test_relative_nested_includes_and_comments(self):
        self.write("models/a.dtd", "<!ENTITY % b SYSTEM '../common/b.dtd'>\n%b;\n<!ELEMENT A (B)>")
        self.write("common/b.dtd", '<!ELEMENT B EMPTY>\n<!ATTLIST B label CDATA "a>b">')
        root = self.write("root.dtd", '<!-- %ignored; -->\n<!ENTITY % a SYSTEM "models/a.dtd">\n%a;')
        bundler = bundle.Bundler()
        result = bundler.expand(root)
        self.assertIn('<!ELEMENT A (B)>', result)
        self.assertIn('<!ELEMENT B EMPTY>', result)
        self.assertIn('label CDATA "a>b"', result)
        self.assertIn('<!-- Begin a.dtd -->', result)
        self.assertNotIn('<!ENTITY', result)
        self.assertEqual(len(bundler.dependencies), 3)

    def test_missing_file_and_cycle(self):
        with self.assertRaisesRegex(bundle.BundleError, r'root.dtd:.*[\s\S]*Cannot read DTD'):
            self.expand('<!ENTITY % a SYSTEM "missing.dtd">%a;')
        self.write('a.dtd', '%a;')
        with self.assertRaisesRegex(bundle.BundleError, 'Include cycle'):
            self.expand('<!ENTITY % a SYSTEM "a.dtd">%a;')

    def test_duplicate_element_reports_both_sources(self):
        self.write('a.dtd', '<!ELEMENT A EMPTY>')
        with self.assertRaisesRegex(bundle.BundleError, 'root.dtd:.*duplicate element A; first at .*a.dtd:'):
            self.expand('<!ENTITY % a SYSTEM "a.dtd">%a;\n<!ELEMENT A ANY>')

    def test_unsupported_and_undefined_entities(self):
        for text in [
            '<!ENTITY % a "<!ELEMENT A EMPTY>">%a;',
            '<!ENTITY % a PUBLIC "public" "a.dtd">%a;',
            '<!ENTITY % a SYSTEM "https://example.com/a.dtd">%a;',
            '<!ENTITY general "text">',
            '<!ELEMENT A (%part;)>',
            '<![INCLUDE[<!ELEMENT A EMPTY>]]>',
            '%undefined;',
            '<!ELEMENT A EMPTY',
            '<!ENTITY % a SYSTEM "a.dtd"><!ENTITY % a SYSTEM "b.dtd">',
        ]:
            with self.subTest(text=text), self.assertRaises(bundle.BundleError):
                self.expand(text)

    def test_cli_failure_preserves_output(self):
        root = self.write('root.dtd', '%undefined;')
        output = self.write('result.dtd', 'previous valid output')
        result = subprocess.run([sys.executable, str(SOURCE / 'config/scripts/bundle_dtd.py'),
                                 str(root), str(output)], capture_output=True, text=True)
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('undefined parameter entity', result.stderr)
        self.assertEqual(output.read_text(), 'previous valid output')

    def test_build_tracks_nested_includes(self):
        self.write('root.dtd', '<!ENTITY % a SYSTEM "models/a.dtd">%a;')
        fragment = self.write('models/a.dtd', '<!ELEMENT A EMPTY>')
        self.write('CMakeLists.txt', '''cmake_minimum_required(VERSION 3.12)
project(DTDDependencyTest NONE)
include("%s/config/cmake/ElffModelDTD.cmake")
elff_bundle_dtd("${CMAKE_CURRENT_SOURCE_DIR}/root.dtd" "${CMAKE_CURRENT_BINARY_DIR}/result.dtd")
''' % SOURCE.as_posix())
        build = self.root / 'build'
        def run(*args):
            result = subprocess.run(['cmake', *map(str, args)], capture_output=True, text=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        run('-S', self.root, '-B', build)
        output = build / 'result.dtd'
        before = output.stat().st_mtime_ns
        run('--build', build)
        self.assertEqual(output.stat().st_mtime_ns, before)
        time.sleep(1.1)  # Also works with build tools using whole-second timestamps.
        # Introduce a new nested dependency through an existing tracked include.
        fragment.write_text('<!ENTITY % b SYSTEM "b.dtd">%b;\n<!ELEMENT A (B)>')
        nested = self.write('models/b.dtd', '<!ELEMENT B EMPTY>')
        run('--build', build)
        self.assertIn('<!ELEMENT B EMPTY>', output.read_text())
        time.sleep(1.1)
        nested.write_text('<!ELEMENT B ANY>')
        run('--build', build)
        self.assertIn('<!ELEMENT B ANY>', output.read_text())


if __name__ == '__main__':
    unittest.main()
