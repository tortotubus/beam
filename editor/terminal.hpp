#pragma once
#include "projection.hpp"
#include <optional>
#include <string>
#include <vector>

struct screen;

namespace editor {
/**
 * @brief Manage terminal input, styled drawing and simple interactive prompts.
 */
class Terminal {
public:
  enum class Mode { Document, Menu, Input };

  enum Key {
    Escape = 27,
    Up = 256,
    Down,
    Left,
    Right,
    Home,
    End,
    Delete,
    BackTab,
    Resize
  };
  /**
   * @brief Enter raw terminal mode and initialize screen and signal handling.
   */
  Terminal();
  /**
   * @brief Restore terminal settings, screen state and previous signal
   * handlers.
   */
  ~Terminal();
  Terminal(const Terminal &) = delete;
  Terminal &operator=(const Terminal &) = delete;
  /**
   * @brief Read a character or decode a special key, resize event or
   * termination request.
   */
  int key();
  /**
   * @brief Refresh the terminal dimensions when the operating system provides
   * them.
   */
  void size();
  /**
   * @brief Draw a clipped viewport with selection, syntax colors, trailing
   * context and status.
   */
  void draw(const std::vector<std::string> &lines, std::size_t selected,
            const std::string &status,
            std::optional<std::size_t> column = std::nullopt,
            std::size_t length = 1,
            std::optional<std::size_t> contextEnd = std::nullopt,
            const std::vector<SyntaxSpan> *syntax = nullptr,
            Mode mode = Mode::Document,
            std::optional<Kind> context = std::nullopt);
  /**
   * @brief Edit a single value, returning its contents on Enter or no value on
   * cancellation.
   */
  std::optional<std::string> input(const std::string &label, std::string value);
  /**
   * @brief Display choices and return the selected index, or no index on
   * cancellation.
   */
  std::optional<std::size_t> menu(const std::string &title,
                                  const std::vector<std::string> &options);
  /**
   * @brief Wrap a diagnostic message into lines that can be browsed in a menu.
   */
  void notice(const std::string &title, const std::string &message);
  int rows = 24, columns = 80;

private:
  /**
   * @brief Draw a status strip and wrapped, mode-specific keyboard hints.
   */
  void footer(const std::string &status, Mode mode,
              std::optional<Kind> context = std::nullopt);

  screen *session = nullptr;
  bool colors = false;
};
} // namespace editor
