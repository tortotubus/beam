#include "editor.hpp"
#include <algorithm>

namespace editor {

void Editor::structural(const Location &location, int command) {
  bool insert = command == '\n';
  if ((insert && location.kind != Kind::Gap) ||
      (!insert && location.kind != Kind::Element)) {
    doc.message =
        insert ? "Select an insertion gap first." : "Select an element first.";
    return;
  }
  if (!insert && (!location.node->parent ||
                  location.node->parent->type != XML_ELEMENT_NODE)) {
    doc.message = "The document root cannot be removed or replaced.";
    return;
  }
  if (command == 'd') {
    doc.removeElement(location.node);
    return;
  }
  auto options = insert ? doc.insertableElements(location.node, location.before)
                        : doc.replacementElements(location.node);
  if (options.empty()) {
    doc.message = "No constructible element is permitted at this position.";
    return;
  }
  auto choice = terminal.menu(
      insert ? "Insert element" : "Replace element (subtree will be discarded)",
      options);
  if (!choice)
    return;
  if (insert)
    doc.insertElement(location.node, location.before, options[*choice]);
  else
    doc.replaceElement(location.node, options[*choice]);
}
void Editor::attribute(xmlNodePtr node, const std::string &name) {
  auto declarations = Constraints(doc.xml()).attributes(node->name);
  if (name.empty()) {
    declarations.erase(std::remove_if(declarations.begin(), declarations.end(),
                                      [&](auto declaration) {
                                        return attributeNode(
                                                   node, declaration->name) !=
                                               nullptr;
                                      }),
                       declarations.end());
  }
  xmlAttributePtr declaration = nullptr;
  if (!name.empty()) {
    for (auto d : declarations)
      if (name == reinterpret_cast<const char *>(d->name))
        declaration = d;
  } else {
    std::vector<std::string> options;
    for (auto d : declarations)
      options.emplace_back(reinterpret_cast<const char *>(d->name));
    if (options.empty()) {
      doc.message = "There are no missing attributes to add.";
      return;
    }
    auto choice = terminal.menu("Add attribute", options);
    if (!choice)
      return;
    declaration = declarations[*choice];
  }
  if (!declaration) {
    doc.message = "Attribute has no DTD declaration.";
    return;
  }
  if (declaration->def == XML_ATTRIBUTE_FIXED) {
    if (!attributeNode(node, declaration->name) && declaration->defaultValue)
      doc.setAttribute(
          node, reinterpret_cast<const char *>(declaration->name),
          reinterpret_cast<const char *>(declaration->defaultValue));
    else
      doc.message = "This attribute is fixed by the DTD.";
    return;
  }
  std::optional<std::string> value;
  if (declaration->tree) {
    std::vector<std::string> options;
    for (auto e = declaration->tree; e; e = e->next)
      options.emplace_back(reinterpret_cast<const char *>(e->name));
    auto choice = terminal.menu("Select value", options);
    if (choice)
      value = options[*choice];
  } else {
    auto current = attributeNode(node, declaration->name);
    value = terminal.input(
        reinterpret_cast<const char *>(declaration->name),
        current ? content(reinterpret_cast<xmlNodePtr>(current)) : "");
  }
  if (value)
    doc.setAttribute(node, reinterpret_cast<const char *>(declaration->name),
                     *value);
}
void Editor::text(const Location &location) {
  auto parent =
      location.kind == Kind::Text ? location.node->parent : location.node;
  if (!Constraints(doc.xml()).textAllowed(parent)) {
    doc.message = "This element does not allow character data.";
    return;
  }
  auto value = terminal.input(
      "Text", location.kind == Kind::Text ? content(location.node) : "");
  if (!value)
    return;
  if (location.kind == Kind::Text)
    doc.setText(location.node, *value);
  else
    doc.insertText(location.node, *value, location.before);
}
void Editor::run() {
  if (doc.message.empty())
    doc.message = "Select an element, attribute, text field, or blank "
                  "insertion gap. ? for help.";
  for (;;) {
    Projection projection(doc);
    cursor.clamp(projection);
    if (projection.locations.empty())
      return;
    auto location = projection.locations[cursor.index()];
    auto contextEnd = projection.lines.size() - 1;
    if (cursor.index() + 1 < projection.locations.size()) {
      auto nextLine = projection.locations[cursor.index() + 1].line;
      contextEnd = nextLine > location.line ? nextLine - 1 : location.line;
    }
    terminal.draw(projection.lines, location.line,
                  (doc.modified() ? "* " : "") + doc.message, location.column,
                  location.length, contextEnd, &projection.styles,
                  Terminal::Mode::Document, location.kind);
    int k = terminal.key();
    std::optional<Cursor::Motion> motion;
    switch (k) {
    case Terminal::Up:
    case 'k':
      motion = Cursor::Motion::Up;
      break;
    case Terminal::Down:
    case 'j':
      motion = Cursor::Motion::Down;
      break;
    case Terminal::Left:
      motion = Cursor::Motion::Left;
      break;
    case Terminal::Right:
      motion = Cursor::Motion::Right;
      break;
    case Terminal::BackTab:
      motion = Cursor::Motion::Previous;
      break;
    case '\t':
      motion = Cursor::Motion::Next;
      break;
    case Terminal::Home:
      motion = Cursor::Motion::First;
      break;
    case Terminal::End:
      motion = Cursor::Motion::Last;
      break;
    }
    if (motion) {
      cursor.move(*motion, projection);
      continue;
    }
    if (k != Terminal::Resize)
      cursor.resetColumn();
    if ((k == 'r' || k == 'd') && location.kind == Kind::Element)
      structural(location, k);
    else if ((k == 'd' || k == 'x') && location.kind == Kind::Attribute) {
      if (k == 'd')
        doc.removeAttribute(location.node, location.attribute);
      else
        doc.resetAttribute(location.node, location.attribute);
    } else if (k == '\r' || k == '\n') {
      if (location.kind == Kind::AttributeGap)
        attribute(location.node);
      else if (location.kind == Kind::Attribute)
        attribute(location.node, location.attribute);
      else if (location.kind == Kind::Text || location.kind == Kind::TextGap)
        text(location);
      else if (location.kind == Kind::Gap)
        structural(location, '\n');
    } else if (k == 's' || k == 19)
      doc.save();
    else if (k == 'v') {
      doc.validate();
      terminal.notice("Validation result", doc.message);
    } else if (k == '!') {
      terminal.notice("Last operation", doc.message);
    } else if (k == '?') {
      terminal.menu(
          "Help",
          {"Up / Down / j / k: previous or next editable line, retaining the "
           "target column",
           "Left / Right: previous or next structural stop on the same line",
           "Tab / Shift-Tab: next or previous structural location across lines",
           "Home / End: first or last location; trailing blank cells allow "
           "insertion",
           "Enter on child gap: insert a permitted element with required "
           "children",
           "Enter on blank inside tag: add a missing "
           "attribute",
           "d / r on element: delete / replace its complete subtree",
           "Enter on attribute: edit its value",
           "d on attribute: remove unless required; x: remove optional or "
           "restore default",
           "Enter on text or its blank placeholder: edit text",
           "Values: Left/Right, Home/End, Backspace/Delete, Ctrl-U clears",
           "Enter commits values; Escape cancels; enum attributes use a menu",
           "s / Ctrl-S: validate and save; v: validate; !: full last message; "
           "q: quit",
           "Required placeholders block saving. Non-ASCII bytes display as "
           "\\xHH.",
           "Long lines are clipped. Value editing scrolls horizontally. No "
           "undo."});
    } else if (k == 'q' || k == 3) {
      if (!doc.modified())
        return;
      auto choice =
          terminal.menu("Unsaved changes", {"Continue editing", "Save and quit",
                                            "Discard changes and quit"});
      if (choice && (*choice == 2 || (*choice == 1 && doc.save())))
        return;
    }
  }
}
} // namespace editor
