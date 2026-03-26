 

#pragma once

#include "elff/config/config.hpp"
#include <iostream>
#include <streambuf>
#include <string>

namespace ELFF
{
class PrefixStreamBuf : public std::streambuf
{
protected:
   std::streambuf *m_dest;
   std::string m_prefix;
   bool m_at_line_start;

   bool WritePrefix();

   int overflow(int ch) override;
   std::streamsize xsputn(const char *s, std::streamsize count) override;
   int sync() override;

public:
   PrefixStreamBuf() : m_dest(nullptr), m_prefix(), m_at_line_start(true) { }

   void SetDestination(std::streambuf *dest) { m_dest = dest; }
   void SetPrefix(const std::string &prefix) { m_prefix = prefix; }
   void ClearPrefix() { m_prefix.clear(); }
   const std::string &GetPrefix() const { return m_prefix; }
   void ResetLineState() { m_at_line_start = true; }
};

class OutStream : public std::ostream
{
protected:
   std::streambuf *m_rdbuf;
   std::ostream *m_tie;
   PrefixStreamBuf m_prefix_buf;

   void Init();

public:
   OutStream(std::ostream &os)
      : std::ostream(NULL), m_rdbuf(nullptr), m_tie(nullptr), m_prefix_buf()
   {
      SetStream(os);
   }
   void SetStream(std::ostream &os)
   {
      m_rdbuf = os.rdbuf();
      m_tie = os.tie();
      m_prefix_buf.SetDestination(m_rdbuf);
      m_prefix_buf.ResetLineState();
      rdbuf(&m_prefix_buf);
      tie(m_tie);
      Init();
   }

   void Enable()
   {
      if (!IsEnabled())
      {
         m_prefix_buf.SetDestination(m_rdbuf);
         m_prefix_buf.ResetLineState();
         rdbuf(&m_prefix_buf);
         tie(m_tie);
      }
   }
   void Disable()
   {
      if (IsEnabled()) { rdbuf(NULL); tie(NULL); }
   }
   bool IsEnabled() const { return (rdbuf() != NULL); }

   void SetPrefix(const std::string &prefix) { m_prefix_buf.SetPrefix(prefix); }
   void ClearPrefix() { m_prefix_buf.ClearPrefix(); }
   const std::string &GetPrefix() const { return m_prefix_buf.GetPrefix(); }
};

extern ELFF_EXPORT OutStream out;
extern ELFF_EXPORT OutStream err;

void SetOutPrefix(const std::string &prefix);
void ClearOutPrefix();
void SetErrPrefix(const std::string &prefix);
void ClearErrPrefix();

std::string MakeParFilename(const std::string &prefix, const int myid,
                            const std::string suffix = "", const int width = 6);
const char* GetEnv(const char* name);

} // namespace ELFF
