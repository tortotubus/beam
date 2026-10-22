#include "elff/io/C/xmlCaseParameters.h"
#include "elff/config/_config.hpp"

#include <cmath>
#include <functional>
#include <locale>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

struct elff_xml_parameters_handle {
  std::map<std::string, double> parameters;
  std::map<std::string, double> coefficients;
};

namespace {
constexpr size_t max_text = 4096;
constexpr size_t max_entries = 4096;
constexpr size_t max_depth = 128;

struct Failure {
  elff_xml_status status;
  size_t line;
  std::string message;
};

[[noreturn]] void invalid(const std::string &message) {
  throw Failure{ELFF_XML_INVALID_VALUE, 0, message};
}

bool letter(char c) {
  return (c >= 'a' && c <= 'z') || (c >= 'A' && c <= 'Z') || c == '_';
}
bool digit(char c) { return c >= '0' && c <= '9'; }
bool space(char c) { return c == ' ' || c == '\t' || c == '\r' || c == '\n'; }
bool identifier(const std::string &s) {
  if (s.empty() || !letter(s[0]))
    return false;
  for (char c : s)
    if (!letter(c) && !digit(c))
      return false;
  return s != "pi";
}
double finite(double value) {
  if (!std::isfinite(value))
    invalid("Nonfinite result or invalid real-valued operation.");
  return value;
}

// Recursive-descent parser; unary -> power ensures -2^2 == -(2^2).
class Expression {
  const std::string &text;
  std::function<double(const std::string &)> resolve;
  size_t pos = 0, depth = 0;
  bool qualified;
  struct Level {
    size_t &depth;
    explicit Level(size_t &d) : depth(d) {
      if (depth >= max_depth)
        invalid("Expression nesting limit exceeded.");
      ++depth;
    }
    ~Level() { --depth; }
  };
  void whitespace() {
    while (pos < text.size() && space(text[pos]))
      ++pos;
  }
  bool take(char c) {
    whitespace();
    if (pos < text.size() && text[pos] == c) {
      ++pos;
      return true;
    }
    return false;
  }
  double number() {
    whitespace();
    const size_t start = pos;
    size_t digits = 0;
    while (pos < text.size() && digit(text[pos])) {
      ++pos;
      ++digits;
    }
    if (pos < text.size() && text[pos] == '.') {
      ++pos;
      while (pos < text.size() && digit(text[pos])) {
        ++pos;
        ++digits;
      }
    }
    if (!digits)
      invalid("Expected a decimal number.");
    if (pos < text.size() && (text[pos] == 'e' || text[pos] == 'E')) {
      ++pos;
      if (pos < text.size() && (text[pos] == '+' || text[pos] == '-'))
        ++pos;
      const size_t exponent = pos;
      while (pos < text.size() && digit(text[pos]))
        ++pos;
      if (exponent == pos)
        invalid("Missing scientific-notation exponent.");
    }
    std::istringstream input(text.substr(start, pos - start));
    input.imbue(std::locale::classic());
    double value;
    if (!(input >> value))
      invalid("Numeric literal is out of range.");
    return finite(value);
  }
  double primary() {
    if (take('(')) {
      double v = sum();
      if (!take(')'))
        invalid("Expected closing parenthesis.");
      return v;
    }
    whitespace();
    if (pos < text.size() && letter(text[pos])) {
      size_t start = pos++;
      while (pos < text.size() && (letter(text[pos]) || digit(text[pos]) || (qualified && text[pos] == '.')))
        ++pos;
      std::string name = text.substr(start, pos - start);
      if (name == "sqrt" && take('(')) {
        double v = sum();
        if (!take(')'))
          invalid("Expected closing parenthesis after sqrt.");
        return finite(std::sqrt(v));
      }
      return name == "pi" ? std::acos(-1.0) : resolve(name);
    }
    return number();
  }
  double power() {
    double v = primary();
    if (take('^'))
      v = finite(std::pow(v, unary()));
    return v;
  }
  double unary() {
    Level level(depth);
    if (take('+'))
      return unary();
    if (take('-'))
      return -unary();
    return power();
  }
  double product() {
    double v = unary();
    for (;;) {
      if (take('*'))
        v = finite(v * unary());
      else if (take('/')) {
        double rhs = unary();
        if (rhs == 0)
          invalid("Division by zero.");
        v = finite(v / rhs);
      } else
        return v;
    }
  }
  double sum() {
    double v = product();
    for (;;) {
      if (take('+'))
        v = finite(v + product());
      else if (take('-'))
        v = finite(v - product());
      else
        return v;
    }
  }

public:
  Expression(const std::string &s, std::function<double(const std::string &)> r,
             bool allow_qualified = false)
      : text(s), resolve(std::move(r)), qualified(allow_qualified) {}
  double evaluate(bool literal = false) {
    double v;
    if (literal) {
      double sign = take('-') ? -1 : 1;
      if (sign == 1)
        take('+');
      v = sign * number();
    } else
      v = sum();
    whitespace();
    if (pos != text.size())
      invalid("Unexpected trailing expression text.");
    return finite(v);
  }
};

#ifdef ELFF_USE_LIBXML2
void check(elff_xml_status status, const elff_xml_error &error) {
  if (status != ELFF_XML_OK)
    throw Failure{status, error.line, error.message};
}
std::string attribute(elff_xml_node_t node, const char *name, bool &present) {
  elff_xml_error error{};
  size_t capacity = 0;
  auto status =
      elff_xml_node_attribute_string(node, name, nullptr, 0, &capacity, &error);
  present = status != ELFF_XML_NOT_FOUND;
  if (!present)
    return {};
  if (status != ELFF_XML_BUFFER_TOO_SMALL)
    check(status, error);
  if (capacity > max_text + 1)
    invalid("Attribute length limit exceeded.");
  std::string value(capacity, '\0');
  check(elff_xml_node_attribute_string(node, name, &value[0], capacity, nullptr,
                                       &error),
        error);
  value.resize(capacity - 1);
  return value;
}
struct Definition {
  std::string text;
  size_t line;
  bool literal;
  int state = 0;
  double value = 0;
};
using Definitions = std::map<std::string, Definition>;
Definitions read_section(elff_xml_node_t root, const char *section_name,
                         const char *entry_name, bool coefficient) {
  elff_xml_error error{};
  elff_xml_node_t section = nullptr;
  auto status =
      elff_xml_node_unique_child(root, section_name, &section, &error);
  if (status == ELFF_XML_NOT_FOUND)
    return {};
  check(status, error);
  size_t count;
  check(elff_xml_node_element_child_count(section, &count, &error), error);
  if (count > max_entries)
    invalid(std::string(section_name) + ": entry limit exceeded.");
  Definitions definitions;
  for (size_t i = 0; i < count; ++i) {
    elff_xml_node_t node;
    check(elff_xml_node_element_child_at(section, i, &node, &error), error);
    size_t line = 0;
    check(elff_xml_node_line(node, &line, &error), error);
    std::string name;
    try {
      const char *tag;
      check(elff_xml_node_name(node, &tag, &error), error);
      if (std::string(tag) != entry_name)
        invalid("Unexpected element in section.");
      bool has_name, has_value, has_expression;
      name = attribute(node, "name", has_name);
      if (!has_name || !identifier(name))
        invalid("Expected a valid name; pi is reserved.");
      std::string value = attribute(node, "value", has_value);
      std::string expression = attribute(node, "expression", has_expression);
      if (coefficient ? (!has_expression || has_value)
                   : (has_value == has_expression))
        invalid(coefficient
                    ? "Coefficient requires expression and no value."
                    : "Parameter requires exactly one of value or expression.");
      if (!definitions
               .emplace(name, Definition{has_value ? value : expression, line,
                                         has_value})
               .second)
        invalid("Duplicate name.");
      size_t children;
      check(elff_xml_node_element_child_count(node, &children, &error), error);
      if (children)
        invalid("Definitions cannot contain child elements.");
    } catch (Failure &failure) {
      failure.line = line;
      failure.message =
          std::string(entry_name) + " '" + name + "': " + failure.message;
      throw;
    }
  }
  return definitions;
}

void resolve_all(Definitions &parameters, Definitions &coefficients,
                 elff_xml_parameters_handle &context) {
  std::function<double(const std::string &, size_t)> resolve;
  resolve = [&](const std::string &name, size_t depth) -> double {
    auto found = parameters.find(name);
    if (found == parameters.end())
      invalid("Undefined parameter '" + name + "'.");
    auto &def = found->second;
    if (def.state == 2)
      return def.value;
    try {
      if (def.state == 1)
        invalid("Cyclic parameter reference to '" + name + "'.");
      if (depth >= max_depth)
        invalid("Parameter dependency limit exceeded.");
      def.state = 1;
      def.value = Expression(def.text, [&](const std::string &dependency) {
                    return resolve(dependency, depth + 1);
                  }).evaluate(def.literal);
      def.state = 2;
      return def.value;
    } catch (Failure &failure) {
      if (!failure.line) {
        failure.line = def.line;
        failure.message = "Parameter '" + name + "': " + failure.message;
      }
      throw;
    }
  };
  for (const auto &entry : parameters)
    context.parameters.emplace(entry.first, resolve(entry.first, 0));
  for (const auto &entry : coefficients) {
    try {
      context.coefficients.emplace(
          entry.first,
          Expression(entry.second.text, [&](const std::string &name) {
            return resolve(name, 0);
          }).evaluate());
    } catch (Failure &failure) {
      failure.line = entry.second.line;
      failure.message = "Coefficient '" + entry.first + "': " + failure.message;
      throw;
    }
  }
}
#endif

elff_xml_status query(elff_xml_parameters_t context, const char *name,
                      double *value, elff_xml_error *error, bool coefficient) {
  elff_xml_error_clear(error);
  if (!context || !name || !*name || !value) {
    elff_xml_error_set(error, 0, 0,
                       "Context, nonempty name and output value are required.");
    return ELFF_XML_INVALID_ARGUMENT;
  }
  try {
    const auto &values =
        coefficient ? context->coefficients : context->parameters;
    auto found = values.find(name);
    if (found == values.end()) {
      elff_xml_error_set(error, 0, 0, "Requested name was not found.");
      return ELFF_XML_NOT_FOUND;
    }
    *value = found->second;
    return ELFF_XML_OK;
  } catch (...) {
    elff_xml_error_set(error, 0, 0, "Unable to query XML parameters.");
    return ELFF_XML_INTERNAL_ERROR;
  }
}
} // namespace

extern "C" elff_xml_status
elff_xml_parameters_create(elff_xml_node_t case_node,
                           elff_xml_parameters_t *out_context,
                           elff_xml_error *error) {
  elff_xml_error_clear(error);
  if (out_context)
    *out_context = nullptr;
  if (!case_node || !out_context) {
    elff_xml_error_set(error, 0, 0,
                       "Case node and output context are required.");
    return ELFF_XML_INVALID_ARGUMENT;
  }
#ifndef ELFF_USE_LIBXML2
  elff_xml_error_set(error, 0, 0, "XML support is unavailable in this build.");
  return ELFF_XML_UNAVAILABLE;
#else
  try {
    elff_xml_error diagnostic{};
    const char *name;
    check(elff_xml_node_name(case_node, &name, &diagnostic), diagnostic);
    if (std::string(name) != "Case")
      invalid("Expected a Case element.");
    auto parameters = read_section(case_node, "Parameters", "Parameter", false);
    auto coefficients =
        read_section(case_node, "Nondimensionalization", "Coefficient", true);
    std::unique_ptr<elff_xml_parameters_handle> context(
        new elff_xml_parameters_handle);
    resolve_all(parameters, coefficients, *context);
    *out_context = context.release();
    return ELFF_XML_OK;
  } catch (const Failure &failure) {
    elff_xml_error_set(error, failure.line, 0, failure.message.c_str());
    return failure.status;
  } catch (...) {
    elff_xml_error_set(error, 0, 0,
                       "Unable to construct XML parameter context.");
    return ELFF_XML_INTERNAL_ERROR;
  }
#endif
}

extern "C" elff_xml_status
elff_xml_parameters_get(elff_xml_parameters_t context, const char *name,
                        double *value, elff_xml_error *error) {
  return query(context, name, value, error, false);
}

extern "C" elff_xml_status
elff_xml_parameters_get_coefficient(elff_xml_parameters_t context,
                                    const char *name, double *value,
                                    elff_xml_error *error) {
  return query(context, name, value, error, true);
}

extern "C" void elff_xml_parameters_destroy(elff_xml_parameters_t context) {
  delete context;
}

extern "C" elff_xml_status
elff_xml_parameters_evaluate(elff_xml_parameters_t context, const char* expression,
                             double* value, elff_xml_error* error)
{
  elff_xml_error_clear(error);
  if (!context || !expression || !value) {
    elff_xml_error_set(error, 0, 0, "Context, expression and output are required.");
    return ELFF_XML_INVALID_ARGUMENT;
  }
  try {
    std::string text(expression);
    if (text.size() > max_text) invalid("Expression length limit exceeded.");
    double result = Expression(text, [&](const std::string& name) {
      const bool coefficient = name.compare(0, 13, "coefficients.") == 0;
      const auto& values = coefficient ? context->coefficients : context->parameters;
      std::string key = coefficient ? name.substr(13) :
        (name.compare(0, 11, "parameters.") == 0 ? name.substr(11) : name);
      auto found = values.find(key);
      if (found == values.end()) invalid("Undefined reference '" + name + "'.");
      return found->second;
    }, true).evaluate();
    *value = result;
    return ELFF_XML_OK;
  } catch (const Failure& failure) {
    elff_xml_error_set(error, failure.line, 0, failure.message.c_str());
    return failure.status;
  } catch (...) {
    elff_xml_error_set(error, 0, 0, "Unable to evaluate render expression.");
    return ELFF_XML_INTERNAL_ERROR;
  }
}
