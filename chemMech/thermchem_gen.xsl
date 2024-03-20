<?xml version="1.0"?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                version="1.0">
  <xsl:output method="text" indent="no" encoding="utf-8"/>

  <xsl:template match="/">
    <xsl:text>= = = = = = = = = = = = = = = = = = = =&#10;</xsl:text>
    <xsl:text>- - - - - SPECIES THERMO DATA - - - - - &#10;</xsl:text>
    <xsl:text>= = = = = = = = = = = = = = = = = = = =&#10;</xsl:text>
    <xsl:apply-templates select="/ctml/speciesData/species/thermo"/>
    <xsl:apply-templates select="/ctml/speciesData/species/thermo/NASA/floatArray"/>
    <xsl:text>= = = = = = = = = = = = = = = = = = = =&#10;</xsl:text>
    <xsl:text>- - - - - REACTION STEP DATA - - - - - &#10;</xsl:text>
    <xsl:text>= = = = = = = = = = = = = = = = = = = =&#10;</xsl:text>
    <xsl:apply-templates select="/ctml/reactionData/reaction/rateCoeff"/>
    <!-- <xsl:text>= = = = = = = = = = = = = = = = = = = =&#10;</xsl:text>
    <xsl:text>- - - - - THREE BODY EFFICIENCY - - - - &#10;</xsl:text>
    <xsl:text>= = = = = = = = = = = = = = = = = = = =&#10;</xsl:text>
    <xsl:apply-templates select="/ctml/reactionData/reaction/rateCoeff/efficiencies"/>
    <xsl:text>= = = = = = = = = = = = = = = = = = = =&#10;</xsl:text>
    <xsl:text>- - - - - FALLOFF REACTIONS - - - - - - &#10;</xsl:text>
    <xsl:text>= = = = = = = = = = = = = = = = = = = =&#10;</xsl:text>
    <xsl:apply-templates select="/ctml/reactionData/reaction/rateCoeff/falloff"/> -->
  </xsl:template>

  <xsl:template match="thermo">
    <xsl:value-of select="concat(../@name, ', 1.0')"/><xsl:text>
</xsl:text>
  </xsl:template>

  <xsl:template match="floatArray">
    <xsl:value-of select="concat(../../../@name, ', ', ../@P0, ', ',
      ../@Tmin, ', ', ../@Tmax, ', ' ,../floatArray/@size)"/>
    <xsl:text>&#10;</xsl:text>
    <xsl:value-of select="normalize-space(.)"/>
    <xsl:text>&#10;</xsl:text>
  </xsl:template>

  <xsl:template match="rateCoeff">
    <xsl:choose>

      <xsl:when test="not(../@type)">
        <xsl:value-of select="concat(../@id, ', ', 'elementary', ', ', ../@reversible)"/>
        <xsl:text>&#10;</xsl:text>
        <xsl:value-of select="../equation"/>
        <xsl:text>&#10;</xsl:text>
        <xsl:value-of select="concat(Arrhenius[1]/A, ', ', Arrhenius[1]/b,
          ', ', Arrhenius[1]/E, ', ', Arrhenius[1]/E/@units)"/>
        <xsl:text>&#10;</xsl:text>
      </xsl:when>

      <xsl:when test="../@type='threeBody'">
        <xsl:value-of select="concat(../@id, ', ', ../@type, ', ', ../@reversible)"/>
        <xsl:text>&#10;</xsl:text>
        <xsl:value-of select="../equation"/>
        <xsl:text>&#10;</xsl:text>
        <xsl:value-of select="concat(Arrhenius[1]/A, ', ', Arrhenius[1]/b,
          ', ', Arrhenius[1]/E, ', ', Arrhenius[1]/E/@units)"/>
        <xsl:text>&#10;</xsl:text>
        <xsl:value-of select="efficiencies"/>
        <xsl:text>&#10;</xsl:text>
      </xsl:when>

      <xsl:otherwise>
        <xsl:value-of select="concat(../@id, ', ', ../@type, ', ', ../@reversible)"/>
        <xsl:text>&#10;</xsl:text>
        <xsl:value-of select="../equation"/>
        <xsl:text>&#10;</xsl:text>
        <xsl:value-of select="concat(Arrhenius[1]/A, ', ', Arrhenius[1]/b,
          ', ', Arrhenius[1]/E, ', ', Arrhenius[1]/E/@units)"/>
        <xsl:text>&#10;</xsl:text>
        <xsl:if test="efficiencies">
          <xsl:value-of select="efficiencies"/>
          <xsl:text>&#10;</xsl:text>
        </xsl:if>

        <xsl:choose>
          <xsl:when test="falloff/@type='Lindemann'">
            <xsl:value-of select="falloff/@type"/>
            <xsl:text>&#10;</xsl:text>
            <xsl:value-of select="concat(Arrhenius[2]/A, ', ', Arrhenius[2]/b, ', ',
              Arrhenius[2]/E, ', ', Arrhenius[2]/E/@units)"/>
            <xsl:text>&#10;</xsl:text>
          </xsl:when>
          <xsl:otherwise>
            <xsl:value-of select="falloff/@type"/>
            <xsl:text>&#10;</xsl:text>
            <xsl:value-of select="concat(Arrhenius[2]/A, ', ', Arrhenius[2]/b, ', ',
              Arrhenius[2]/E, ', ', Arrhenius[2]/E/@units)"/>
            <xsl:text>&#10;</xsl:text>
            <xsl:value-of select="falloff"/>
            <xsl:text>&#10;</xsl:text>
          </xsl:otherwise>
        </xsl:choose>

      </xsl:otherwise>
    </xsl:choose>
    <xsl:text>***************************************&#10;</xsl:text>
  </xsl:template>

  <!-- <xsl:template match="efficiencies">
    <xsl:value-of select="concat(../../@id, ', ', ../../@type)"/>
    <xsl:text>&#10;</xsl:text>
    <xsl:value-of select="../../equation"/>
    <xsl:text>&#10;</xsl:text>
    <xsl:value-of select="../efficiencies"/>
    <xsl:text>&#10;</xsl:text>
    <xsl:text>***************************************&#10;</xsl:text>
  </xsl:template> -->

  <!-- <xsl:template match="falloff">
    <xsl:value-of select="concat(../../@id, ', ', @type)"/>
    <xsl:text>&#10;</xsl:text>
    <xsl:value-of select="../../equation"/>
    <xsl:text>&#10;</xsl:text>
    <xsl:value-of select="concat(../Arrhenius[2]/A, ', ', ../Arrhenius[2]/b, ', ',
      ../Arrhenius[2]/E, ', ', ../Arrhenius[2]/E/@units)"/>
      <xsl:text>&#10;</xsl:text>
    <xsl:choose>
      <xsl:when test="@type='Lindemann'">
      </xsl:when>
      <xsl:otherwise test="@type='Troe'">
        <xsl:value-of select="../falloff"/>
        <xsl:text>&#10;</xsl:text>
      </xsl:otherwise>
  </xsl:choose>
    <xsl:text>***************************************&#10;</xsl:text>
  </xsl:template> -->

</xsl:stylesheet>
