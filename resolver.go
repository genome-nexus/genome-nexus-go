package genome_nexus_annotator_go

import (
	"regexp"
	"strconv"
	"strings"

	gnapi "github.com/genome-nexus/genome-nexus-go-api-client/genome-nexus-public-api"
	tempotype "github.mskcc.org/cdsi/cdsi-protobuf/tempo/generated/v2/go"
)

/*
Utility functions for resolving response (variant annotation) from genome nexus
against CVR-specific data (tempotype.Event)
*/

const (
	defaultBuildNumber string = "37"
	defaultStrand      string = "+"
)

var (
	dbeventRsidRegex      = regexp.MustCompile("^(rs\\d*)$")
	proteinPositionRegex  = regexp.MustCompile("p.[A-Za-z]([0-9]*).*$")
	validNucleotidesRegex = regexp.MustCompile("^([ATGC]*)$")
)

func resolveSiftScore(canonicalTranscript gnapi.TranscriptConsequenceSummary) string {
	if canonicalTranscript != (gnapi.TranscriptConsequenceSummary{}) && canonicalTranscript.SiftScore != nil {
		if *canonicalTranscript.SiftScore == float64(int64(*canonicalTranscript.SiftScore)) {
			return strconv.FormatFloat(*canonicalTranscript.SiftScore, 'f', 1, 64)
		}
		return strconv.FormatFloat(*canonicalTranscript.SiftScore, 'g', 3, 64)
	}
	return ""
}

func resolveSiftPrediction(canonicalTranscript gnapi.TranscriptConsequenceSummary) string {
	if canonicalTranscript != (gnapi.TranscriptConsequenceSummary{}) && canonicalTranscript.SiftPrediction != nil {
		return *canonicalTranscript.SiftPrediction
	}
	return ""
}

func resolvePolyphenScore(canonicalTranscript gnapi.TranscriptConsequenceSummary) string {
	if canonicalTranscript != (gnapi.TranscriptConsequenceSummary{}) && canonicalTranscript.PolyphenScore != nil {
		if *canonicalTranscript.PolyphenScore == float64(int64(*canonicalTranscript.PolyphenScore)) {
			return strconv.FormatFloat(*canonicalTranscript.PolyphenScore, 'f', 1, 64)
		}
		return strconv.FormatFloat(*canonicalTranscript.PolyphenScore, 'g', 3, 64)
	}
	return ""
}

func resolvePolyphenPrediction(canonicalTranscript gnapi.TranscriptConsequenceSummary) string {
	if canonicalTranscript != (gnapi.TranscriptConsequenceSummary{}) && canonicalTranscript.PolyphenPrediction != nil {
		return *canonicalTranscript.PolyphenPrediction
	}
	return ""
}

func getCanonicalTranscript(gnResponse gnapi.VariantAnnotation) gnapi.TranscriptConsequenceSummary {
	// TODO: Handle case where we return empty struct
	if gnResponse.AnnotationSummary != nil &&
		gnResponse.AnnotationSummary.TranscriptConsequences != nil &&
		len(gnResponse.AnnotationSummary.TranscriptConsequences) > 0 {
		return gnResponse.AnnotationSummary.TranscriptConsequences[0]
	}
	return gnapi.TranscriptConsequenceSummary{}
}

func resolveHugoSymbol(canonicalTranscript gnapi.TranscriptConsequenceSummary) string {
	if canonicalTranscript != (gnapi.TranscriptConsequenceSummary{}) && canonicalTranscript.HugoGeneSymbol != nil {
		return *canonicalTranscript.HugoGeneSymbol
	}
	return "FAILED"
}

func resolveEntrezGeneId(canonicalTranscript gnapi.TranscriptConsequenceSummary) string {
	// TODO is this okay to return 0 if canonicalTranscript.EntrezGeneId is empty?
	if canonicalTranscript == (gnapi.TranscriptConsequenceSummary{}) || canonicalTranscript.EntrezGeneId == nil {
		return ""
	}
	return *canonicalTranscript.EntrezGeneId
}

func resolveAssemblyName(gnResponse gnapi.VariantAnnotation) string {
	if gnResponse.AssemblyName == nil {
		return defaultBuildNumber
	}
	return *gnResponse.AssemblyName
}

func resolveChromosome(gnResponse gnapi.VariantAnnotation, genomicLocation gnapi.GenomicLocation) string {
	if gnResponse.AnnotationSummary != nil {
		chromosome, ok := gnResponse.AnnotationSummary.GenomicLocation.GetChromosomeOk()
		if ok {
			return *chromosome
		}
	}
	return genomicLocation.Chromosome
}

func resolveStart(gnResponse gnapi.VariantAnnotation, genomicLocation gnapi.GenomicLocation) string {
	// Try to get start from VariantAnnotation
	if gnResponse.AnnotationSummary != nil {
		start, ok := gnResponse.AnnotationSummary.GenomicLocation.GetStartOk()
		if ok {
			return strconv.Itoa(int(*start))
		}
	}

	// Else try to get start from GenomicLocation
	start, ok := genomicLocation.GetStartOk()
	if ok {
		return strconv.Itoa(int(*start))
	}

	// Return empty string otherwise
	return ""
}

func resolveEnd(gnResponse gnapi.VariantAnnotation, genomicLocation gnapi.GenomicLocation) string {
	// Try to get end from VariantAnnotation
	if gnResponse.AnnotationSummary != nil {
		end, ok := gnResponse.AnnotationSummary.GenomicLocation.GetEndOk()
		if ok {
			return strconv.Itoa(int(*end))
		}
	}

	// Else try to get end from GenomicLocation
	end, ok := genomicLocation.GetEndOk()
	if ok {
		return strconv.Itoa(int(*end))
	}

	// Return empty string otherwise
	return ""
}

func resolveStrandSign(gnResponse gnapi.VariantAnnotation) string {
	if gnResponse.AnnotationSummary != nil && gnResponse.AnnotationSummary.StrandSign != nil {
		return *gnResponse.AnnotationSummary.StrandSign
	}
	return defaultStrand
}

func resolveVariantClassification(canonicalTranscript gnapi.TranscriptConsequenceSummary, event tempotype.Event) string {
	if canonicalTranscript != (gnapi.TranscriptConsequenceSummary{}) &&
		canonicalTranscript.VariantClassification != nil {
		return *canonicalTranscript.VariantClassification
	}
	return event.VariantClassification
}

// NOTE: Why don't we get the variant type from resolveCVRVariantType if AnnotationSummary.VariantType is empty?
func resolveVariantType(gnResponse gnapi.VariantAnnotation) string {
	if gnResponse.AnnotationSummary != nil && gnResponse.AnnotationSummary.VariantType != nil {
		return *gnResponse.AnnotationSummary.VariantType
	}
	return ""
}

func resolveDbSnpRs(gnResponse gnapi.VariantAnnotation) string {
	if gnResponse.ColocatedVariants != nil && len(gnResponse.ColocatedVariants) > 0 {
		for _, cv := range gnResponse.ColocatedVariants {
			// TODO check if DbSnpId is nil?
			var match = dbeventRsidRegex.FindStringSubmatch(*cv.DbSnpId)
			if match != nil {
				return match[0]
			}
		}
	}
	return ""
}

func resolveHgvsc(canonicalTranscript gnapi.TranscriptConsequenceSummary) string {
	if canonicalTranscript != (gnapi.TranscriptConsequenceSummary{}) && canonicalTranscript.Hgvsc != nil {
		return *canonicalTranscript.Hgvsc
	}
	return ""
}

func resolveHgvsp(canonicalTranscript gnapi.TranscriptConsequenceSummary) string {
	if canonicalTranscript != (gnapi.TranscriptConsequenceSummary{}) && canonicalTranscript.Hgvsp != nil {
		return *canonicalTranscript.Hgvsp
	}
	return ""
}

func resolveHgvspShort(canonicalTranscript gnapi.TranscriptConsequenceSummary) string {
	if canonicalTranscript != (gnapi.TranscriptConsequenceSummary{}) && canonicalTranscript.HgvspShort != nil {
		return *canonicalTranscript.HgvspShort
	}
	return ""
}

func resolveTranscriptId(canonicalTranscript gnapi.TranscriptConsequenceSummary) string {
	if canonicalTranscript != (gnapi.TranscriptConsequenceSummary{}) {
		transcriptId, ok := canonicalTranscript.GetTranscriptIdOk()
		if ok {
			return *transcriptId
		}
	}
	return ""
}

func resolveRefSeq(canonicalTranscript gnapi.TranscriptConsequenceSummary) string {
	if canonicalTranscript != (gnapi.TranscriptConsequenceSummary{}) && canonicalTranscript.RefSeq != nil {
		return *canonicalTranscript.RefSeq
	}
	return ""
}

func resolveProteinPosStart(canonicalTranscript gnapi.TranscriptConsequenceSummary) string {
	if canonicalTranscript != (gnapi.TranscriptConsequenceSummary{}) {
		start, ok := canonicalTranscript.ProteinPosition.GetStartOk()
		if ok {
			return strconv.Itoa(int(*start))
		}
	}
	return ""
}

func resolveProteinPosEnd(canonicalTranscript gnapi.TranscriptConsequenceSummary) string {
	if canonicalTranscript != (gnapi.TranscriptConsequenceSummary{}) {
		end, ok := canonicalTranscript.ProteinPosition.GetEndOk()
		if ok {
			return strconv.Itoa(int(*end))
		}
	}
	return ""
}

func resolveCodonChange(canonicalTranscript gnapi.TranscriptConsequenceSummary) string {
	if canonicalTranscript != (gnapi.TranscriptConsequenceSummary{}) && canonicalTranscript.CodonChange != nil {
		return *canonicalTranscript.CodonChange
	}
	return ""
}

func resolveConsequence(canonicalTranscript gnapi.TranscriptConsequenceSummary) string {
	if canonicalTranscript != (gnapi.TranscriptConsequenceSummary{}) && canonicalTranscript.ConsequenceTerms != nil {
		return *canonicalTranscript.ConsequenceTerms
	}
	return ""
}

func resolveProteinPosition(canonicalTranscript gnapi.TranscriptConsequenceSummary) string {
	if canonicalTranscript != (gnapi.TranscriptConsequenceSummary{}) {
		var proteinPosStart string = resolveProteinPosStart(canonicalTranscript)
		var hgvspShort string = resolveHgvspShort(canonicalTranscript)
		if proteinPosStart != "" {
			return proteinPosStart
		} else if hgvspShort != "" {
			// try extracting from hgvspShort if proteinPosStart empty
			var match = proteinPositionRegex.FindStringSubmatch(hgvspShort)
			if match != nil && match[1] != "" {
				return match[1]
			}
		}
	}
	// NOTE: In the original method, it tries to get proteinPosition from the mutationRecord's additionalProperties
	return ""
}

func resolveExon(canonicalTranscript gnapi.TranscriptConsequenceSummary) string {
	if canonicalTranscript != (gnapi.TranscriptConsequenceSummary{}) && canonicalTranscript.Exon != nil {
		return *canonicalTranscript.Exon
	}
	return ""
}

func resolveReferenceAllele(gnResponse gnapi.VariantAnnotation, event tempotype.Event) string {
	if gnResponse.AnnotationSummary != nil {
		referenceAllele, ok := gnResponse.AnnotationSummary.GenomicLocation.GetReferenceAlleleOk()
		if ok {
			return *referenceAllele
		}
	}
	return event.ReferenceAllele
}

func resolveTumorSeqAllele(gnResponse gnapi.VariantAnnotation, event tempotype.Event) string {
	if gnResponse.AnnotationSummary != nil {
		variantAllele, ok := gnResponse.AnnotationSummary.GenomicLocation.GetVariantAlleleOk()
		if ok {
			return *variantAllele
		}
	}
	return resolveTumorSeqAlleleFromInput(event.ReferenceAllele, event.TumorSeqAllele1, event.TumorSeqAllele2)
}

func resolveTumorSeqAlleleFromInput(referenceAllele string, tumorSeqAllele1 string, tumorSeqAllele2 string) string {
	// Note below functionality comes from here: https://github.com/cBioPortal/cbioportal/blob/2c107f919e5ef959379b94e058acf41341571202/maf/src/main/java/org/mskcc/cbio/maf/MafUtil.java#L739-L774
	// Sanity check tumor seq allele 1 and 2 for valid/non-null values
	if (tumorSeqAllele1 == "" || strings.EqualFold(tumorSeqAllele1, "NA")) &&
		(tumorSeqAllele2 == "" || strings.EqualFold(tumorSeqAllele2, "NA")) {
		// Cannot resolve this case
		return ""
	}

	// Resolve tumor seq allele
	if tumorSeqAllele1 == "" || strings.EqualFold(tumorSeqAllele1, "NA") || tumorSeqAllele1 == referenceAllele {
		return tumorSeqAllele2
	} else if tumorSeqAllele2 == "" || strings.EqualFold(tumorSeqAllele2, "NA") || tumorSeqAllele2 == referenceAllele {
		return tumorSeqAllele1
	} else if variantContainsAmbiguousTumorSeqAllele(referenceAllele, tumorSeqAllele1, tumorSeqAllele2) {
		if tumorSeqAllele2 == "-" {
			return tumorSeqAllele1
		}
		return tumorSeqAllele2
	}
	return tumorSeqAllele1
}

// This function returns true when one tumorSeqAllele is a '-' and one contains a valid nucleotide pattern
// In all other cases it returns false, including if both tumorSeqAlleles contain valid nucleotide patterns
func variantContainsAmbiguousTumorSeqAllele(referenceAllele string, tumorSeqAllele1 string, tumorSeqAllele2 string) bool {
	// Tumor seq allele 1 or 2 is null type or equal to ref allele - return false
	if tumorSeqAllele1 == "" || strings.EqualFold(tumorSeqAllele1, "NA") || tumorSeqAllele1 == referenceAllele ||
		tumorSeqAllele2 == "" || strings.EqualFold(tumorSeqAllele2, "NA") || tumorSeqAllele2 == referenceAllele {
		return false
	}

	// Returns true if one of the tumorSeqAlleles is a '-' and one contains a valid nucleotide pattern
	return (tumorSeqAllele1 == "-" || tumorSeqAllele2 == "-") &&
		(validNucleotidesRegex.MatchString(strings.ToUpper(tumorSeqAllele1)) ||
			validNucleotidesRegex.MatchString(strings.ToUpper(tumorSeqAllele2)))
}

func resolveRefAndTumorSeqAlleles(gnResponse gnapi.VariantAnnotation, event tempotype.Event, stripMatchingBases string) (string, string, string) {
	resolvedReferenceAllele := resolveReferenceAllele(gnResponse, event)
	resolvedTumorSeqAllele1 := event.ReferenceAllele
	resolvedTumorSeqAllele2 := resolveTumorSeqAllele(gnResponse, event)

	// Get tumorSeqAllele from input, it could be from tumorSeqAllele2 or tumorSeqAllele1
	// Logic is here: https://github.com/cBioPortal/cbioportal/blob/master/core/src/main/java/org/mskcc/cbio/maf/MafUtil.java#L811
	resolvedTumorSeqAlleleFromInput := resolveTumorSeqAlleleFromInput(event.ReferenceAllele, event.TumorSeqAllele1, event.TumorSeqAllele2)

	// If keep all allele bases
	if stripMatchingBases == "none" {
		// If keep all allele bases, referenceAllele and tumorSeqAllele1 would be the input value
		// TumorSeqAllele2 would be the resolved result that was sent to genome nexus server.
		resolvedReferenceAllele = event.ReferenceAllele
		resolvedTumorSeqAllele2 = resolvedTumorSeqAlleleFromInput
	} else if stripMatchingBases == "first" {
		// If strip first allele bases, first check if referenceAllele has matching bases, and length > 1 to ensure stripping is valid.
		// Then remove the first matching allele base
		if len(event.ReferenceAllele) > 1 && event.ReferenceAllele != resolvedReferenceAllele {
			resolvedReferenceAllele = event.ReferenceAllele[1:len(event.ReferenceAllele)]
		} else {
			resolvedReferenceAllele = event.ReferenceAllele
		}
		// If resolvedTumorSeqAlleleFromInput equals resolvedTumorSeqAllele2, it means there are matching allele bases in the input and needs to do stripping
		// Check if length > 1 to ensure stripping is valid
		if resolvedTumorSeqAlleleFromInput != resolvedTumorSeqAllele2 && len(resolvedTumorSeqAlleleFromInput) > 1 {
			resolvedTumorSeqAllele2 = resolvedTumorSeqAlleleFromInput[1:len(event.ReferenceAllele)]
		}
	}

	// NOTE: In the original genome-nexus-annotation-pipeline implementation, we check whether
	// tumorSeqAllele1 equals to referenceAllele or tumorSeqAllele2 and, if yes, assign its value to tumorSeqAllele1.
	// Otherwise tumorSeqAllele1 stays as is.
	// We are no longer performing this check because tumorSeqAllele1 is initialized with the value of RefAllele
	// and is not modified through reading of the MAF or any other operation by the time this code is reached.
	// Meaning, tumorSeqAllele1 is always equal to the value of RefAllele.
	resolvedTumorSeqAllele1 = resolvedReferenceAllele

	return resolvedReferenceAllele, resolvedTumorSeqAllele1, resolvedTumorSeqAllele2
}
