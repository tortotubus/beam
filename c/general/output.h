#pragma once

#ifdef __cplusplus
extern "C"
{
#endif

/**
 * @brief Set a line prefix for the ELFF standard output stream.
 *
 * Passing `NULL` clears the prefix.
 *
 * @param prefix Prefix string to prepend to each new output line
 */
void elff_set_out_prefix(const char *prefix);

/**
 * @brief Clear any line prefix on the ELFF standard output stream.
 */
void elff_clear_out_prefix(void);

/**
 * @brief Set a line prefix for the ELFF standard error stream.
 *
 * Passing `NULL` clears the prefix.
 *
 * @param prefix Prefix string to prepend to each new error line
 */
void elff_set_err_prefix(const char *prefix);

/**
 * @brief Clear any line prefix on the ELFF standard error stream.
 */
void elff_clear_err_prefix(void);

#ifdef __cplusplus
}
#endif
