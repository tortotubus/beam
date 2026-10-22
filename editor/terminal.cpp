#include "terminal.hpp"
#include <algorithm>
#include <csignal>
#include <cstdio>
#include <curses.h>
#include <poll.h>
#include <stdexcept>
#include <unistd.h>

namespace editor {
namespace {
volatile sig_atomic_t stopped = 0;
/**
 * @brief Record a signal for orderly handling by the input loop.
 */
void stop(int) { stopped = 1; }
const int signals[] = {SIGINT, SIGTERM, SIGHUP, SIGQUIT, SIGTSTP};
struct sigaction previous[5];
/**
 * @brief Locate the preceding UTF-8 code-point boundary.
 */
std::size_t previousChar(const std::string &text, std::size_t pos) {
  if (pos)
    --pos;
  while (pos && (static_cast<unsigned char>(text[pos]) & 0xc0) == 0x80)
    --pos;
  return pos;
}
/**
 * @brief Locate the following UTF-8 code-point boundary.
 */
std::size_t nextChar(const std::string &text, std::size_t pos) {
  if (pos < text.size())
    ++pos;
  while (pos < text.size() &&
         (static_cast<unsigned char>(text[pos]) & 0xc0) == 0x80)
    ++pos;
  return pos;
}
} // namespace
Terminal::Terminal() {
  if (!isatty(STDIN_FILENO) || !isatty(STDOUT_FILENO))
    throw std::runtime_error("Interactive mode requires a terminal.");
  stopped = 0;
  for (int i = 0; i < 5; ++i)
    sigaction(signals[i], nullptr, &previous[i]);
  session = newterm(nullptr, stdout, stdin);
  if (!session)
    throw std::runtime_error("Cannot initialize ncurses; check TERM.");
  raw();
  noecho();
  keypad(stdscr, TRUE);
  wtimeout(stdscr, 100);
  set_escdelay(40);
  curs_set(0);
  colors = has_colors() && start_color() == OK;
  if (colors) {
    const short background = use_default_colors() == OK ? -1 : COLOR_BLACK;
    init_pair(1, COLOR_CYAN, background);
    init_pair(2, COLOR_GREEN, background);
    init_pair(3, COLOR_YELLOW, background);
    init_pair(4, COLOR_WHITE, COLOR_RED);
    init_pair(5, COLOR_WHITE, COLOR_BLUE);
    init_pair(6, COLOR_BLACK, COLOR_CYAN);
  }
  struct sigaction action{};
  action.sa_handler = stop;
  sigemptyset(&action.sa_mask);
  for (int i = 0; i < 5; ++i)
    sigaction(signals[i], &action, nullptr);
  size();
}
Terminal::~Terminal() {
  endwin();
  delscreen(session);
  for (int i = 0; i < 5; ++i)
    sigaction(signals[i], &previous[i], nullptr);
}
void Terminal::size() { getmaxyx(stdscr, rows, columns); }
int Terminal::key() {
  for (;;) {
    if (stopped)
      throw std::runtime_error("Terminated; unsaved changes were discarded.");
    int k = wgetch(stdscr);
    switch (k) {
    case KEY_UP:
      return Up;
    case KEY_DOWN:
      return Down;
    case KEY_LEFT:
      return Left;
    case KEY_RIGHT:
      return Right;
    case KEY_HOME:
      return Home;
    case KEY_END:
      return End;
    case KEY_DC:
      return Delete;
    case KEY_BTAB:
      return BackTab;
    case KEY_BACKSPACE:
      return 127;
    case KEY_ENTER:
      return '\n';
    case KEY_RESIZE:
      size();
      return Resize;
    case ERR: {
      pollfd pfd{STDIN_FILENO, POLLIN, 0};
      poll(&pfd, 1, 0);
      if (pfd.revents & (POLLHUP | POLLERR | POLLNVAL))
        throw std::runtime_error("Terminal disconnected.");
      break;
    }
    default:
      if (k >= 0 && k < 256)
        return k;
      // Unknown function keys are not editable characters.
      break;
    }
  }
}
void Terminal::draw(const std::vector<std::string> &lines, std::size_t selected,
                    const std::string &status,
                    std::optional<std::size_t> column, std::size_t length,
                    std::optional<std::size_t> contextEnd,
                    const std::vector<SyntaxSpan> *syntax, Mode mode,
                    std::optional<Kind> context) {
  size();
  curs_set(0);
  erase();
  int height = std::max(1, rows - 3);
  // Include read-only lines after the cursor, without scrolling the cursor
  // out of view. Closing tags have no editing stops of their own.
  auto last =
      contextEnd
          ? std::max(selected, std::min(*contextEnd, selected + height - 1))
          : selected;
  std::size_t first =
      last >= static_cast<std::size_t>(height) ? last - height + 1 : 0;
  for (int row = 0; row < height; ++row) {
    auto i = first + row;
    if (i < lines.size()) {
      auto width = static_cast<std::size_t>(std::max(1, columns - 1));
      auto offset =
          i == selected && column && *column >= width ? *column - width + 1 : 0;
      auto shown = lines[i].substr(std::min(offset, lines[i].size()), width);
      // Apply styles after clipping to preserve projection coordinates.
      std::vector<chtype> styles(shown.size(), A_NORMAL);
      if (syntax) {
        for (const auto &span : *syntax) {
          if (span.line != i)
            continue;
          auto begin = std::max(offset, span.column);
          auto end = std::min(offset + shown.size(), span.column + span.length);
          for (auto position = begin; position < end; ++position)
            styles[position - offset] =
                colors ? COLOR_PAIR(static_cast<int>(span.color)) : A_NORMAL;
        }
      }
      if (i == selected) {
        auto begin = column ? *column - offset : 0;
        auto end =
            column ? std::min(shown.size(), begin + length) : shown.size();
        // Selection overrides syntax and placeholder colors.
        for (auto position = begin; position < end; ++position)
          styles[position] = A_REVERSE;
      }
      for (std::size_t begin = 0; begin < shown.size();) {
        auto end = begin + 1;
        while (end < shown.size() && styles[end] == styles[begin])
          ++end;
        attrset(styles[begin]);
        mvaddnstr(row, static_cast<int>(begin), shown.data() + begin,
                  static_cast<int>(end - begin));
        begin = end;
      }
    }
  }
  footer(status, mode, context);
  refresh();
}
void Terminal::footer(const std::string &status, Mode mode,
                      std::optional<Kind> context) {
  const int top = std::max(0, rows - 3);
  const int width = std::max(1, columns - 1);
  for (int row = top; row < rows; ++row) {
    attrset(A_NORMAL);
    move(row, 0);
    clrtoeol();
  }
  attrset(colors ? COLOR_PAIR(5) | A_BOLD : A_REVERSE);
  mvhline(top, 0, ' ', width);
  mvaddnstr(top, 0, visible(status).c_str(), width);

  using Hint = std::pair<const char *, const char *>;
  std::vector<Hint> hints;
  switch (mode) {
  case Mode::Document:
    if (context == Kind::Element)
      hints = {{"d", "Delete"}, {"r", "Replace"}};
    else if (context == Kind::AttributeGap)
      hints = {{"Enter", "Add attribute"}};
    else if (context == Kind::Gap)
      hints = {{"Enter", "Insert element"}};
    else {
      hints = {{"Enter", "Edit"}};
      if (context == Kind::Attribute) {
        hints.push_back({"d", "Delete"});
        hints.push_back({"x", "Reset"});
      }
    }
    for (const Hint hint : {Hint{"s", "Save"},
                            {"?", "Help"},
                            {"q", "Quit"},
                            {"Arrows", "Move"},
                            {"Tab", "Next"},
                            {"v", "Validate"},
                            {"!", "Details"}})
      hints.push_back(hint);
    break;
  case Mode::Menu:
    hints = {{"Enter", "Select"},
             {"Esc", "Cancel"},
             {"Up/Down", "Move"},
             {"Tab", "Next"}};
    break;
  case Mode::Input:
    hints = {{"Enter", "Apply"},        {"Esc", "Cancel"},
             {"Ctrl-U", "Clear"},       {"Left/Right", "Move"},
             {"Home/End", "Start/end"}, {"Bksp/Del", "Delete"}};
    break;
  }
  int row = top + 1, column = 0;
  for (const auto &hint : hints) {
    const std::string key = " " + std::string(hint.first) + " ";
    const std::string label = " " + std::string(hint.second) + "  ";
    const int length = static_cast<int>(key.size() + label.size());
    if (column && column + length > width) {
      ++row;
      column = 0;
    }
    if (row >= rows)
      break;
    // Keep hints whole when the terminal is too narrow to show every action.
    if (length > width)
      continue;
    attrset(colors ? COLOR_PAIR(6) | A_BOLD : A_REVERSE | A_BOLD);
    mvaddnstr(row, column, key.c_str(), static_cast<int>(key.size()));
    attrset(A_NORMAL);
    addnstr(label.c_str(), static_cast<int>(label.size()));
    column += length;
  }
  attrset(A_NORMAL);
}
std::optional<std::string> Terminal::input(const std::string &label,
                                           std::string value) {
  std::size_t cursor = value.size();
  for (;;) {
    size();
    auto prefix = visible(label) + ": ";
    auto shown = prefix + visible(value);
    auto column = prefix.size() + visible(value.substr(0, cursor)).size();
    auto width = static_cast<std::size_t>(std::max(1, columns - 1));
    auto offset = column >= width ? column - width + 1 : 0;
    footer("", Mode::Input);
    const int inputRow = std::max(0, rows - 3);
    attrset(colors ? COLOR_PAIR(5) | A_BOLD : A_REVERSE);
    move(inputRow, 0);
    addnstr(shown.c_str() + offset, static_cast<int>(width));
    move(inputRow, static_cast<int>(column - offset));
    attrset(A_NORMAL);
    curs_set(1);
    refresh();
    int k = key();
    if (k == Escape || k == 3)
      return std::nullopt;
    if (k == '\r' || k == '\n')
      return value;
    if (k == Left)
      cursor = previousChar(value, cursor);
    else if (k == Right)
      cursor = nextChar(value, cursor);
    else if (k == Home || k == 1)
      cursor = 0;
    else if (k == End || k == 5)
      cursor = value.size();
    else if ((k == 127 || k == 8) && cursor) {
      auto start = previousChar(value, cursor);
      value.erase(start, cursor - start);
      cursor = start;
    } else if (k == Delete && cursor < value.size())
      value.erase(cursor, nextChar(value, cursor) - cursor);
    else if (k == 21) {
      value.clear();
      cursor = 0;
    } else if (k >= 32 && k < 256 && k != 127)
      value.insert(cursor++, 1, static_cast<char>(k));
  }
}
void Terminal::notice(const std::string &title, const std::string &message) {
  size();
  std::vector<std::string> parts;
  auto width = static_cast<std::size_t>(std::max(1, columns - 1));
  std::size_t start = 0;
  while (start < message.size()) {
    auto end = message.find('\n', start);
    if (end == std::string::npos)
      end = message.size();
    auto line = visible(message.substr(start, end - start));
    for (std::size_t pos = 0; pos < line.size(); pos += width)
      parts.push_back(line.substr(pos, width));
    start = end + 1;
  }
  if (parts.empty())
    parts.push_back("No diagnostic.");
  menu(title, parts);
}
std::optional<std::size_t>
Terminal::menu(const std::string &title,
               const std::vector<std::string> &options) {
  if (options.empty())
    return std::nullopt;
  std::vector<std::string> lines;
  for (const auto &item : options)
    lines.push_back(visible(item));
  std::size_t selected = 0;
  for (;;) {
    draw(lines, selected,
         title + "  [" + std::to_string(selected + 1) + "/" +
             std::to_string(options.size()) + "]",
         std::nullopt, 1, std::nullopt, nullptr, Mode::Menu);
    int k = key();
    if (k == Escape || k == 3 || k == 'q')
      return std::nullopt;
    if (k == '\r' || k == '\n')
      return selected;
    if ((k == Up || k == 'k' || k == BackTab) && selected)
      --selected;
    if ((k == Down || k == 'j' || k == '\t') && selected + 1 < options.size())
      ++selected;
  }
}
} // namespace editor
