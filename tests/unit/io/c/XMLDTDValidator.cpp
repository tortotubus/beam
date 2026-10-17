#include <gtest/gtest.h>

#include <elff/config/_config.hpp>
#include <elff/io/C/xmlDTDValidator.h>

#include <cstring>
#include <fstream>
#include <iterator>
#include <string>

namespace {

constexpr char valid_elff_xml[] = R"XML(
<ELFF schemaVersion="1" dimensions="2" units="nondimensional">
  <Models>
    <Model id="fibre">
      <EulerBeamHuang>
        <Geometry length="1" nodes="65"/>
        <Material bendingRigidity="0.001" linearDensity="1.5"/>
        <InitialLine anchorEnd="right">
          <Anchor x="0" y="0" z="0"/>
          <Direction x="1" y="0" z="0"/>
        </InitialLine>
        <BoundaryConditions left="free" right="simple"/>
      </EulerBeamHuang>
    </Model>
  </Models>
</ELFF>
)XML";

const std::string& elff_case_dtd()
{
  static const std::string dtd = [] {
    std::ifstream input(ELFF_CASE_DTD_PATH, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(input),
                       std::istreambuf_iterator<char>());
  }();
  return dtd;
}

class ELFFCaseDTDTest : public ::testing::Test
{
protected:
  void SetUp() override
  {
    ASSERT_FALSE(elff_case_dtd().empty())
      << "Could not read " << ELFF_CASE_DTD_PATH;
  }
};

#ifdef ELFF_USE_LIBXML2

elff_xml_dtd_status validate(const char* xml, char* diagnostic, size_t capacity)
{
  const std::string& dtd = elff_case_dtd();
  return elff_xml_dtd_validate_memory(xml,
                                      std::strlen(xml),
                                      dtd.data(),
                                      dtd.size(),
                                      diagnostic,
                                      capacity);
}

elff_xml_dtd_status validate_elff_subtree(const char* xml,
                                          char* diagnostic,
                                          size_t capacity)
{
  const std::string& dtd = elff_case_dtd();
  return elff_xml_dtd_validate_subtree_memory(xml,
                                              std::strlen(xml),
                                              "Case",
                                              "ELFF",
                                              dtd.data(),
                                              dtd.size(),
                                              diagnostic,
                                              capacity);
}

TEST_F(ELFFCaseDTDTest, AcceptsVersionOneElffSection)
{
  char diagnostic[256] = "not cleared";

  EXPECT_EQ(validate(valid_elff_xml, diagnostic, sizeof(diagnostic)),
            ELFF_XML_DTD_VALID)
    << diagnostic;
  EXPECT_STREQ(diagnostic, "");
}

TEST_F(ELFFCaseDTDTest, RejectsModelWithoutId)
{
  constexpr char xml[] = R"XML(
<ELFF schemaVersion="1" dimensions="2" units="nondimensional">
  <Models><Model><EulerBeamHuang/></Model></Models>
</ELFF>
)XML";
  char diagnostic[256] = {};

  EXPECT_EQ(validate(xml, diagnostic, sizeof(diagnostic)),
            ELFF_XML_DTD_INVALID_DOCUMENT);
  EXPECT_NE(diagnostic[0], '\0');
}

TEST_F(ELFFCaseDTDTest, RejectsOtherSchemaVersions)
{
  std::string xml = valid_elff_xml;
  const std::string version = R"(schemaVersion="1")";
  const size_t version_position = xml.find(version);
  char diagnostic[256] = {};

  ASSERT_NE(version_position, std::string::npos);
  xml[version_position + version.size() - 2] = '2';
  EXPECT_EQ(validate(xml.c_str(), diagnostic, sizeof(diagnostic)),
            ELFF_XML_DTD_INVALID_DOCUMENT);
  EXPECT_NE(diagnostic[0], '\0');
}

TEST_F(ELFFCaseDTDTest, RejectsUnknownModelElement)
{
  constexpr char xml[] = R"XML(
<ELFF schemaVersion="1" dimensions="2" units="nondimensional">
  <Models><Model id="fibre"><UnknownBeam/></Model></Models>
</ELFF>
)XML";

  EXPECT_EQ(validate(xml, nullptr, 0), ELFF_XML_DTD_INVALID_DOCUMENT);
}

TEST_F(ELFFCaseDTDTest, DistinguishesMalformedXmlFromDtdViolation)
{
  constexpr char xml[] = "<ELFF><Models></ELFF>";
  char diagnostic[256] = {};

  EXPECT_EQ(validate(xml, diagnostic, sizeof(diagnostic)),
            ELFF_XML_DTD_MALFORMED_XML);
  EXPECT_NE(diagnostic[0], '\0');
}

TEST_F(ELFFCaseDTDTest, ReportsMalformedDtd)
{
  constexpr char malformed_dtd[] = "<!ELEMENT ELFF (Models)";
  char diagnostic[256] = {};

  EXPECT_EQ(elff_xml_dtd_validate_memory(valid_elff_xml,
                                         sizeof(valid_elff_xml) - 1,
                                         malformed_dtd,
                                         sizeof(malformed_dtd) - 1,
                                         diagnostic,
                                         sizeof(diagnostic)),
            ELFF_XML_DTD_MALFORMED_DTD);
  EXPECT_NE(diagnostic[0], '\0');
}

TEST_F(ELFFCaseDTDTest, RejectsNullInput)
{
  const std::string& dtd = elff_case_dtd();
  char diagnostic[64] = {};

  EXPECT_EQ(elff_xml_dtd_validate_memory(nullptr,
                                         0,
                                         dtd.data(),
                                         dtd.size(),
                                         diagnostic,
                                         sizeof(diagnostic)),
            ELFF_XML_DTD_INVALID_ARGUMENT);
  EXPECT_NE(diagnostic[0], '\0');
}

TEST_F(ELFFCaseDTDTest, ValidatesOnlySelectedElffSubtree)
{
  constexpr char xml[] = R"XML(
<Case formatVersion="1">
  <ELFF schemaVersion="1" dimensions="2" units="nondimensional">
    <Models>
      <Model id="fibre">
        <EulerBeamHuang>
          <Geometry length="1" nodes="65"/>
          <Material bendingRigidity="0.001" linearDensity="1.5"/>
          <InitialLine anchorEnd="right">
            <Anchor x="0" y="0"/><Direction x="1" y="0"/>
          </InitialLine>
          <BoundaryConditions left="free" right="simple"/>
        </EulerBeamHuang>
      </Model>
    </Models>
  </ELFF>
  <Basilisk><CaseSpecificElement arbitrary="yes"/></Basilisk>
</Case>
)XML";
  char diagnostic[256] = "not cleared";

  EXPECT_EQ(validate_elff_subtree(xml, diagnostic, sizeof(diagnostic)),
            ELFF_XML_DTD_VALID)
    << diagnostic;
  EXPECT_STREQ(diagnostic, "");
}

TEST_F(ELFFCaseDTDTest, RejectsInvalidSelectedElffSubtree)
{
  constexpr char xml[] = R"XML(
<Case>
  <ELFF schemaVersion="1" dimensions="2" units="nondimensional">
    <Models/>
  </ELFF>
  <Basilisk><Anything/></Basilisk>
</Case>
)XML";
  char diagnostic[256] = {};

  EXPECT_EQ(validate_elff_subtree(xml, diagnostic, sizeof(diagnostic)),
            ELFF_XML_DTD_INVALID_DOCUMENT);
  EXPECT_NE(diagnostic[0], '\0');
}

TEST_F(ELFFCaseDTDTest, ReportsMissingAndDuplicateSelectedSubtrees)
{
  constexpr char missing[] = "<Case><Basilisk/></Case>";
  constexpr char duplicate[] = "<Case><ELFF/><ELFF/></Case>";

  EXPECT_EQ(validate_elff_subtree(missing, nullptr, 0),
            ELFF_XML_DTD_SUBTREE_NOT_FOUND);
  EXPECT_EQ(validate_elff_subtree(duplicate, nullptr, 0),
            ELFF_XML_DTD_SUBTREE_NOT_UNIQUE);
}

TEST_F(ELFFCaseDTDTest, ProvidesStableStatusStrings)
{
  EXPECT_STREQ(elff_xml_dtd_status_string(ELFF_XML_DTD_VALID), "valid");
  EXPECT_STREQ(elff_xml_dtd_status_string(ELFF_XML_DTD_MALFORMED_DTD),
               "malformed DTD");
  EXPECT_STREQ(elff_xml_dtd_status_string(ELFF_XML_DTD_SUBTREE_NOT_FOUND),
               "subtree not found");
  EXPECT_STREQ(elff_xml_dtd_status_string(ELFF_XML_DTD_SUBTREE_NOT_UNIQUE),
               "subtree not unique");
  EXPECT_STREQ(elff_xml_dtd_status_string(
                 static_cast<elff_xml_dtd_status>(-1)),
               "unknown status");
}

#else

TEST_F(ELFFCaseDTDTest, ReportsWhenLibXml2IsDisabled)
{
  const std::string& dtd = elff_case_dtd();

  EXPECT_EQ(elff_xml_dtd_validate_memory(valid_elff_xml,
                                         sizeof(valid_elff_xml) - 1,
                                         dtd.data(),
                                         dtd.size(),
                                         nullptr,
                                         0),
            ELFF_XML_DTD_UNAVAILABLE);
}

#endif

} // namespace
