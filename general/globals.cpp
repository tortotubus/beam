
#include "elff/config/config.hpp"
#include "elff/general/globals.hpp"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdlib>  // getenv

namespace ELFF
{

OutStream out(std::cout);
OutStream err(std::cerr);

namespace Internal
{
bool elff_out_initialized = false;
bool elff_err_initialized = false;
}

bool PrefixStreamBuf::WritePrefix()
{
   if (m_dest == nullptr || m_prefix.empty() || !m_at_line_start)
   {
      return true;
   }

   const std::streamsize written =
      m_dest->sputn(m_prefix.data(), static_cast<std::streamsize>(m_prefix.size()));
   if (written != static_cast<std::streamsize>(m_prefix.size()))
   {
      return false;
   }

   m_at_line_start = false;
   return true;
}

int PrefixStreamBuf::overflow(int ch)
{
   if (ch == traits_type::eof())
   {
      return traits_type::not_eof(ch);
   }

   if (m_dest == nullptr)
   {
      return traits_type::eof();
   }

   if (!WritePrefix())
   {
      return traits_type::eof();
   }

   const int result = m_dest->sputc(static_cast<char>(ch));
   if (result == traits_type::eof())
   {
      return traits_type::eof();
   }

   m_at_line_start = (ch == '\n');
   return result;
}

std::streamsize PrefixStreamBuf::xsputn(const char *s, std::streamsize count)
{
   if (m_dest == nullptr)
   {
      return 0;
   }

   std::streamsize written = 0;
   for (; written < count; ++written)
   {
      if (overflow(static_cast<unsigned char>(s[written])) == traits_type::eof())
      {
         break;
      }
   }
   return written;
}

int PrefixStreamBuf::sync()
{
   return (m_dest == nullptr || m_dest->pubsync() == 0) ? 0 : -1;
}

void OutStream::Init()
{
   if (this == &ELFF::out)
   {
      Internal::elff_out_initialized = true;
   }
   else if (this == &ELFF::err)
   {
      Internal::elff_err_initialized = true;
   }
}

void SetOutPrefix(const std::string &prefix)
{
   ELFF::out.SetPrefix(prefix);
}

void ClearOutPrefix()
{
   ELFF::out.ClearPrefix();
}

void SetErrPrefix(const std::string &prefix)
{
   ELFF::err.SetPrefix(prefix);
}

void ClearErrPrefix()
{
   ELFF::err.ClearPrefix();
}

std::string MakeParFilename(const std::string &prefix, const int myid,
                            const std::string suffix, const int width)
{
   std::stringstream fname;
   fname << prefix << std::setw(width) << std::setfill('0') << myid << suffix;
   return fname.str();
}


const char *GetEnv(const char* name)
{
   return getenv(name);
}

}
