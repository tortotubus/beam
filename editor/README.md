# Structured XML terminal editor

A standalone C++17 editor using libxml2, ncurses, and the C++ standard library.
It has no dependency on the parent ELFF build.

Build and run from the repository root:

```sh
cmake -S editor -B editor/build
cmake --build editor/build
editor/build/xml-editor editor/examples/case.xml
editor/build/xml-editor --check editor/examples/case.xml
```

Install your platform's libxml2 and ncurses development packages, CMake, and a
C++17 compiler, or enter the repository's `nix develop` environment.
Run on a POSIX terminal with a valid `TERM` setting. `--check` needs no terminal and
returns 0 for a valid document, 1 for a parse/validation error, or 2 for bad usage.
The example edits in place; copy it elsewhere to experiment. The executable
contains the repository's `io/DTD/ELFF/ELFF-v2.dtd`, embedded by CMake at build
time. Rebuilding picks up schema changes. No DTD file is needed at runtime.

To start a document, pass a path that does not yet exist:

```sh
editor/build/xml-editor /path/to/new-case.xml
```

The parent directory must already exist. The editor constructs a minimal `Case`
with required children and pending required attributes. Fill those fields and
press `s` to validate and create the file. Until then, the document exists only
in memory; quitting without saving leaves no file. `--check` only checks existing
files and never creates a document. Existing malformed or invalid files are
reported as errors, never replaced with a new skeleton.

## Editing

The bottom status strip shows the last result or the active input field.
Two rows of highlighted key hints follow the selected location and switch
between document commands, menu
selection, and value editing. Hints wrap as the terminal resizes; narrow windows
show the primary actions first. Press `?` in the document for the full key list.
Menus also show the selected item number and total choices.

XML delimiters (`<`, `>`, `</`, `/>`) are cyan, tag names green, and attribute
names yellow. Values retain the terminal's default color. Empty text placeholders
and pending required attributes have a red background. Selection overrides this
with the same reverse-video highlight as other stops. Red text placeholders indicate
an empty editable field, not a DTD validation error: `#PCDATA` may be empty.

The XML tree is authoritative. Tags and XML punctuation are read-only. The
highlight moves between elements, editable text, attributes, and blank insertion
gaps. Attributes appear inline inside the opening tag; selection highlights the
individual value. Insertion stops occupy a trailing blank cell on the preceding line, highlighted
only when selected; they do not add rows. Required-field markers and insertion
cells are UI annotations, never written to the document. A stop appears only
when the insertion menu has at least one permitted, constructible element.
Missing declared attributes (including `#IMPLIED`) have a separate blank stop
inside the opening tag, just before `>` or `/>`. Enter at this gap is the only
way to add a missing attribute; Enter on an existing value edits that value. For example, `<Parameter name="myParameter" />` offers `value`
and `expression` from the ELFF v2 DTD.
Formatting whitespace does not create duplicate stops. Text-only elements
render inline (`<Iterations>2</Iterations>`). Empty elements that permit text
show a selectable blank between paired tags, such as `<LinearDensity> </LinearDensity>`;
this blank is not stored in the XML tree. Move Right onto it and press Enter
to add text. Tag names offer only subtree deletion and replacement; Enter on
a tag name does nothing. Clearing the text restores the blank field. Other childless elements
render as `<tag />` only when their DTD permits neither text nor child elements.
Empty containers use paired tags, such as `<Parameters> </Parameters>`, with
the child insertion stop between the tags. Text and attributes are escaped in
the projection.

| Key | Action |
| --- | --- |
| Up / Down, k / j | Previous / next line with editable stops, nearest to the target column |
| Left / Right | Previous / next structural stop on the same line |
| Tab / Shift-Tab | Next / previous structural location across lines |
| Home / End | First / last location |
| Enter on a value or text placeholder | Edit or add its contents |
| Enter on child gap | Insert a permitted element |
| Enter on blank inside tag | Add a missing declared attribute |
| d / r on element | Delete / replace the entire subtree |
| d on attribute | Delete the explicit attribute unless it is `#REQUIRED` |
| x on attribute | Remove an optional attribute, or restore its DTD default |
| s / Ctrl-S | Validate and save |
| v | Validate and display diagnostics |
| ! | Display the full last operation message |
| ? | Help |
| q / Ctrl-C | Quit; offer save/discard if modified |

Vertical movement remembers the target column across consecutive Up/Down presses
and skips lines containing only read-only syntax. The viewport includes trailing
read-only context, such as closing tags below the last selectable location,
while keeping the selection visible. Left/Right stop at the ends of
the current line; they do not wrap. Use Enter to edit a selected value, then
Left/Right move between characters within that value.

Menus use Up/Down and Enter; Escape cancels. Value editing supports Left/Right,
Home/End, Backspace/Delete, and Ctrl-U to clear. Enter commits and Escape cancels.
Use Enter on the blank inside the opening tag to add optional attributes. Enumerations and NOTATION
attributes present their declared values as a menu. Fixed attribute values cannot be edited, but `d` can remove their explicit
occurrence because `#FIXED` does not require an attribute to be present. Deleting
an attribute with a DTD default removes it from the tree and saved XML; its DTD
default still applies and will be materialized again when the file is reopened.
The attribute insertion menu can restore omitted attributes, using the declared
value for `#FIXED`. All deletions still undergo document validation (for example,
removing an ID that is referenced elsewhere is rejected). Empty CDATA values are valid: explicitly committing an empty required
CDATA attribute resolves its placeholder.

## DTD enforcement

The editor uses one schema: the compiled-in ELFF v2 DTD. Existing input must
have a `Case` root and validate against that schema. A DOCTYPE is optional;
source DOCTYPE declarations are discarded and external DTD paths are not loaded.
Internal declarations cannot override the embedded contract. Defaults are
materialized from the embedded DTD. Saved XML has no DOCTYPE and does not depend
on a DTD file being installed beside it. Other consumers must supply the ELFF v2
schema themselves if they require DTD validation.

XML parsing uses `XML_PARSE_NONET` without external DTD loading. Namespaces and
general entity references remain unsupported. DTD validity does not establish
application-level validity: numeric ranges and relationships described only in
DTD comments are not enforced by this editor.

`constraints.cpp` matches child-element sequences directly against
`xmlElementContent`. It handles sequence, choice, `?`, `*`, `+`, EMPTY, ANY, and
mixed content. Memoized endpoint sets and a finite repetition worklist handle
nullable repetitions without looping. Menus preserve sibling order and show
only substitutions/insertions accepted by the parent's content model.

New elements receive the smallest finite required-child construction found by
DTD relaxation. Choice ties follow declaration order. Optional children are
omitted; required attributes without defaults get explicit placeholders. Models
with no finite construction, or constructions exceeding 4096 elements, are not
offered. Defaults and fixed attributes are populated from the DTD.

Every mutation runs on a copied tree with a freshly parsed copy of the embedded
DTD, preserving the declaration links needed for required-attribute checks. A failed
edit restores the previous tree. Structural constraints always apply. Completed
documents also undergo full libxml2 DTD validation on every edit, including
ID/IDREF, enumeration, and fixed-value checks. While required attribute
placeholders exist, the document is explicitly incomplete: global constraints
such as IDREF resolution can remain unsatisfied. Individual entered attributes
are still validated. Resolving the final placeholder requires full validity;
resolve other reference problems first if necessary. Saving is blocked until
all placeholders are resolved and full DTD validation succeeds.

Save writes UTF-8 to a temporary file in the destination directory and flushes
it. For existing documents, it renames the temporary file over the original.
For new documents, it atomically creates the destination without overwriting a
file that appeared during editing. New files have permissions `0600`. The file's permission bits are preserved; the
inode is replaced (hard links, ownership, ACLs and extended attributes are not
preserved). Symlink destinations are refused. There is no backup or external
change detection. DTD declarations are not editable.

## Scope and implementation

- `document.*`: loading/creation, validated editing operations, private transactions
  and placeholder bookkeeping, and saving.
- `schema.hpp.in`: build-time embedding of the single ELFF v2 schema.
- `constraints.*`: DTD lookup, content-model matching, required-child construction.
- `projection.*`: builds plain display lines, tree-backed cursor locations, and
  color spans without opening a terminal.
- `cursor.*`: moves between projection locations and retains the preferred
  column for vertical navigation; independent of terminal key codes.
- `terminal.*`: ncurses screen lifetime, keys, viewport, value fields and menus.
- `editor.*`: coordinates key bindings, menus, navigation, and document operations.
- `main.cpp`: interactive entry point and `--check`.

Each editor iteration builds a `Projection` from the document, clamps the
`Cursor` to its available locations, and draws the projection. Navigation keys
map to `Cursor::Motion`; editing commands gather values or choices, then call
`Document` operations such as `insertElement`, `setText`, or `removeAttribute`.
Only `Document` mutates the tree and tracks pending attributes. Its transaction
callback and bookkeeping are private.

`Document::xml()` exposes borrowed libxml2 pointers for inspection. Successful
mutations replace the entire tree, invalidating all earlier node pointers and
projections; rejected edits preserve them. The editor rebuilds its projection
before the next input. Insertion indices refer to actual child nodes, including
whitespace and comments, rather than just child elements. The projection and
cursor can be exercised without linking the terminal implementation.

This is intended for small configuration documents. Each edit copies the tree;
menus and rendering inspect DTD models. There is no undo, search, mouse support,
DTD editing, or free-form XML mode. Root deletion/replacement is disabled.
Namespaces and general entity references are rejected explicitly; predefined
XML escapes work normally. Comments and processing instructions are preserved
but read-only. Recursive DTDs with finite expansions are supported.

Whitespace and mixed-content text are preserved on save: indentation is a view,
not a rewrite of text nodes. Existing text newlines remain intact, but value
entry is single-line. Non-ASCII and control bytes are shown as ASCII escapes
(`\xHH`, `\n`, `\t`) to keep terminal control and column calculations predictable.
UTF-8 input is accepted and checked for valid XML 1.0 characters; character
navigation follows UTF-8 code points, not grapheme clusters. Long source lines
are clipped; the selected line scrolls horizontally to reveal the selected
attribute. Value fields also scroll horizontally. Very small terminal windows are
not a supported layout. Terminal settings restore on normal exit, exceptions,
and handled termination signals; SIGKILL cannot be handled. SIGTSTP terminates
rather than suspends, discarding unsaved edits.

The project intentionally inspects libxml2 DTD structures. libxml2 2.15 marks
these interfaces deprecated; their deprecation warnings are suppressed for
this executable. A future libxml2 release that removes them will require an
adaptation.
