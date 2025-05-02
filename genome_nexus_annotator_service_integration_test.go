package genome_nexus_annotator_go

import (
	"encoding/json"
	"os"
	"strings"
	"testing"

	tt "github.mskcc.org/cdsi/cdsi-protobuf/tempo/generated/v1/go"
)

// testdata generated using genome-nexus-annotation-pipeline (Java) annotated MAF
// with the following properties:
//
// -Dgenomenexus.enrichment_fields=annotation_summary,sift,polyphen (application.properties)
// -i<isoform_override=mskcc (cmdline arg)
//
// test can break if genome nexus is updated, in which case rerun (example below) and save new testdata
// $JAVA_BINARY -Dgenomenexus.enrichment_fields=annotation_summary,sift,polyphen -jar annotator.jar -f <inputmaf> -o <outputmaf> -i mskcc

const (
	mutationRecordsJSON = "testdata/tempo_message.annotated.json"
	//mutationRecordsJSON = "testdata/tempo_message.100k.annotated.json"
	isoformOverrideString = "mskcc"
)

type Testset struct {
	Records []tt.TempoMessage `json:"records"`
}

func TestAnnotateMutations(t *testing.T) {
	testset := readTestSetJSON(t, mutationRecordsJSON) // this is a TempoMessage struct with X events for testing (completely populated)
	// for every tempo message in the set set
	for _, tm := range testset.Records {
		// make a deep copy
		testtm := tt.TempoMessage{}
		b, _ := json.Marshal(tm)
		json.Unmarshal(b, &testtm)
		// annotate and compare
		AnnotateTempoMessageEvents(isoformOverrideString, &testtm)
		for i2, event := range(testtm.Events) {
			assertNoError(t, tm.CmoSampleId, event, tm.Events[i2])
		}
	}
}

func readTestSetJSON(t testing.TB, mafFile string) Testset {
	mafData, err := os.ReadFile(mafFile)
	if err != nil {
		t.Fatalf("Failed to read MAF records from JSON %q: %v", mafFile, err)
	}

	//var tempoMessage tt.TempoMessage
	var tempoMessages Testset
	err = json.Unmarshal(mafData, &tempoMessages)
	if err != nil {
		t.Fatalf("Unable to unmarshal test data")
	}
	return tempoMessages
}

var fieldMap = map[int]string{
	0:  "Chromosome",
	1:  "StartPosition",
	2:  "EndPosition",
	3:  "ReferenceAllele",
	4:  "TumorSeqAllele1",
	5:  "TumorSeqAllele2",
	6:  "Strand",
	7:  "NcbiBuild",
	8:  "HugoSymbol",
	9:  "VariantClassification",
	10: "VariantType",
	11: "DbsnpRs",
	12: "Hgvsp",
	13: "HgvspShort",
	14: "Hgvsc",
	15: "TranscriptId",
	16: "Refseq",
	17: "Center",
	18: "Consequence",
	19: "DbsnpValStatus",
	20: "MatchNormSampleBarcode",
	21: "MatchNormSeqAllele1",
	22: "MatchNormSeqAllele2",
	23: "VerificationStatus",
	24: "ValidationStatus",
	25: "MutationStatus",
	26: "SequencingPhase",
	27: "SequencingSource",
	28: "ValidationMethod",
	29: "Score",
	30: "BamFile",
	31: "Sequencer",
	32: "TRefCount",
	33: "TAltCount",
	34: "NRefCount",
	35: "NAltCount",
	36: "ProteinPosion",
	37: "Codons",
	38: "ExonNumber",
	39: "PolyphenPrediction",
	40: "PolyphenScore",
	41: "SiftPrediction",
	42: "SiftScore",
	43: "GenomicLocationExplanation",
	44: "AnnotationStatus",
	45: "EntrezGeneId",
}

func assertNoError(t testing.TB, s string, e *tt.Event, ae *tt.Event) {
	if !strings.EqualFold(e.Chromosome, ae.Chromosome) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[0], e.Hgvsc, ae.Hgvsc)
	}
	if !strings.EqualFold(e.StartPosition, ae.StartPosition) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[1], e.StartPosition, ae.StartPosition)
	}
	if !strings.EqualFold(e.EndPosition, ae.EndPosition) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[2], e.EndPosition, ae.EndPosition)
	}
	if !strings.EqualFold(e.ReferenceAllele, ae.ReferenceAllele) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[3], e.ReferenceAllele, ae.ReferenceAllele)
	}
	if !strings.EqualFold(e.TumorSeqAllele1, ae.TumorSeqAllele1) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[4], e.TumorSeqAllele1, ae.TumorSeqAllele1)
	}
	if !strings.EqualFold(e.TumorSeqAllele2, ae.TumorSeqAllele2) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[5], e.TumorSeqAllele2, ae.TumorSeqAllele2)
	}
	if !strings.EqualFold(e.Strand, ae.Strand) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[6], e.Strand, ae.Strand)
	}
	if !strings.EqualFold(e.NcbiBuild, ae.NcbiBuild) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[7], e.NcbiBuild, ae.NcbiBuild)
	}
	if !strings.EqualFold(e.HugoSymbol, ae.HugoSymbol) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[8], e.HugoSymbol, ae.HugoSymbol)
	}
	if !strings.EqualFold(e.VariantClassification, ae.VariantClassification) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[9], e.VariantClassification, ae.VariantClassification)
	}
	if !strings.EqualFold(e.VariantType, ae.VariantType) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[10], e.VariantType, ae.VariantType)
	}
	if !strings.EqualFold(e.DbsnpRs, ae.DbsnpRs) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[11], e.DbsnpRs, ae.DbsnpRs)
	}
	if !strings.EqualFold(e.Hgvsp, ae.Hgvsp) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[12], e.Hgvsp, ae.Hgvsp)
	}
	if !strings.EqualFold(e.HgvspShort, ae.HgvspShort) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[13], e.HgvspShort, ae.HgvspShort)
	}
	if !strings.EqualFold(e.Hgvsc, ae.Hgvsc) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[14], e.Hgvsc, ae.Hgvsc)
	}
	if !strings.EqualFold(e.TranscriptId, ae.TranscriptId) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[15], e.TranscriptId, ae.TranscriptId)
	}
	if !strings.EqualFold(e.Refseq, ae.Refseq) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[16], e.Refseq, ae.Refseq)
	}
	if !strings.EqualFold(e.Center, ae.Center) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[17], e.Center, ae.Center)
	}
	if !strings.EqualFold(e.Consequence, ae.Consequence) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[18], e.Consequence, ae.Consequence)
	}
	if !strings.EqualFold(e.DbsnpValStatus, ae.DbsnpValStatus) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[19], e.DbsnpValStatus, ae.DbsnpValStatus)
	}
	if !strings.EqualFold(e.MatchedNormSampleBarcode, ae.MatchedNormSampleBarcode) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[20], e.MatchedNormSampleBarcode, ae.MatchedNormSampleBarcode)
	}
	if !strings.EqualFold(e.MatchNormSeqAllele1, ae.MatchNormSeqAllele1) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[21], e.MatchNormSeqAllele1, ae.MatchNormSeqAllele1)
	}
	if !strings.EqualFold(e.MatchNormSeqAllele2, ae.MatchNormSeqAllele2) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[22], e.MatchNormSeqAllele2, ae.MatchNormSeqAllele2)
	}
	if !strings.EqualFold(e.VerificationStatus, ae.VerificationStatus) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[23], e.VerificationStatus, ae.VerificationStatus)
	}
	if !strings.EqualFold(e.ValidationStatus, ae.ValidationStatus) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[24], e.ValidationStatus, ae.ValidationStatus)
	}
	if !strings.EqualFold(e.MutationStatus, ae.MutationStatus) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[25], e.MutationStatus, ae.MutationStatus)
	}
	if !strings.EqualFold(e.SequencingPhase, ae.SequencingPhase) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[26], e.SequencingPhase, ae.SequencingPhase)
	}
	if !strings.EqualFold(e.SequencingSource, ae.SequencingSource) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[27], e.SequencingSource, ae.SequencingSource)
	}
	if !strings.EqualFold(e.ValidationMethod, ae.ValidationMethod) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[28], e.ValidationMethod, ae.ValidationMethod)
	}
	if !strings.EqualFold(e.Score, ae.Score) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[29], e.Score, ae.Score)
	}
	if !strings.EqualFold(e.BamFile, ae.BamFile) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[30], e.BamFile, ae.BamFile)
	}
	if !strings.EqualFold(e.Sequencer, ae.Sequencer) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[31], e.Sequencer, ae.Sequencer)
	}
	if !strings.EqualFold(e.TRefCount, ae.TRefCount) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[32], e.TRefCount, ae.TRefCount)
	}
	if !strings.EqualFold(e.TAltCount, ae.TAltCount) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[33], e.TAltCount, ae.TAltCount)
	}
	if !strings.EqualFold(e.NRefCount, ae.NRefCount) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[34], e.NRefCount, ae.NRefCount)
	}
	if !strings.EqualFold(e.NAltCount, ae.NAltCount) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[35], e.NAltCount, ae.NAltCount)
	}
	if !strings.EqualFold(e.ProteinPosition, ae.ProteinPosition) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[36], e.ProteinPosition, ae.ProteinPosition)
	}
	if !strings.EqualFold(e.Codons, ae.Codons) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[37], e.Codons, ae.Codons)
	}
	if !strings.EqualFold(e.ExonNumber, ae.ExonNumber) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[38], e.ExonNumber, ae.ExonNumber)
	}
	if !strings.EqualFold(e.PolyphenPrediction, ae.PolyphenPrediction) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[39], e.PolyphenPrediction, ae.PolyphenPrediction)
	}
	if !strings.EqualFold(e.PolyphenScore, ae.PolyphenScore) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[40], e.PolyphenScore, ae.PolyphenScore)
	}
	if !strings.EqualFold(e.SiftPrediction, ae.SiftPrediction) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[41], e.SiftPrediction, ae.SiftPrediction)
	}
	if !strings.EqualFold(e.SiftScore, ae.SiftScore) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[42], e.SiftScore, ae.SiftScore)
	}
	if !strings.EqualFold(e.AnnotationStatus, ae.AnnotationStatus) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[44], e.AnnotationStatus, ae.AnnotationStatus)
	}
	if !strings.EqualFold(e.EntrezGeneId, ae.EntrezGeneId) {
		t.Errorf("patient: %q; field: %q; expected %q but got %q", s, fieldMap[45], e.EntrezGeneId, ae.EntrezGeneId)
	}
}
