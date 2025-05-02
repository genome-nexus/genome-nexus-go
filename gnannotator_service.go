package genome_nexus_annotator_go

import (
	"context"
	"fmt"
	"strconv"

	gnapi "github.com/genome-nexus/genome-nexus-go-api-client/genome-nexus-public-api"
	tt "github.mskcc.org/cdsi/cdsi-protobuf/tempo/generated/v1/go"
)

func AnnotateTempoMessageEvents(isoformOverrideSource string, tm *tt.TempoMessage) error {
	// prepare genomic locations for annotation
	genomicLocations := make([]gnapi.GenomicLocation, 0)
	for _, a := range tm.Events {
		loc := getGenomicLocation(a)
		genomicLocations = append(genomicLocations, loc)
	}

	variantAnnotations, err := getVariantAnnotations(isoformOverrideSource, genomicLocations)
	if err != nil {
		return err
	}
	if len(variantAnnotations) != len(genomicLocations) {
		return fmt.Errorf("Returned number of variant annotations does not match genomic locations submitted")
	}

	// construct a list of annotated events using genome nexus response
	for i, variantAnnotation := range variantAnnotations {
		err = mapResponseToEvent(variantAnnotation, genomicLocations[i], tm.Events[i])
		if err != nil {
			return err
		}
	}
	return nil
}

func getGenomicLocation(e *tt.Event) gnapi.GenomicLocation {
	start, _ := strconv.ParseInt(e.StartPosition, 10, 32)
	end, _ := strconv.ParseInt(e.StartPosition, 10, 32)

	gloc := gnapi.NewGenomicLocation(e.Chromosome,
		int32(start),
		int32(end),
		e.ReferenceAllele,
		resolveTumorSeqAlleleFromInput(e.ReferenceAllele, e.TumorSeqAllele1, e.TumorSeqAllele2))
	return *gloc
}

func getVariantAnnotations(isoformOverrideSource string, genomicLocations []gnapi.GenomicLocation) ([]gnapi.VariantAnnotation, error) {
	var token = ""
	fields := make([]string, 0)
	fields = append(fields, "annotation_summary")
	configuration := gnapi.NewConfiguration()
	apiClient := gnapi.NewAPIClient(configuration)

	ctx := context.WithValue(context.Background(), gnapi.ContextServerIndex, 1)
	x := apiClient.AnnotationControllerApi.FetchVariantAnnotationByGenomicLocationPOST(ctx).GenomicLocations(genomicLocations)
	x = x.Fields(fields)
	x = x.IsoformOverrideSource(isoformOverrideSource)
	x = x.Token(token)

	variantAnnotations, r, err := x.Execute()
	if err != nil {
		return variantAnnotations, fmt.Errorf("Error calling Genome Nexus annotation service: %v\nFull HTTP response: %v\n", err, r)
	}
	return variantAnnotations, nil
}

func mapResponseToEvent(variantAnnotation gnapi.VariantAnnotation, genomicLocation gnapi.GenomicLocation, event *tt.Event) error {
	if !*variantAnnotation.SuccessfullyAnnotated {
		event.AnnotationStatus = "FAILED"
		return fmt.Errorf("Unsuccessful variant annotation for genomicLocation: %v", variantAnnotation)
	}
	canonicalTranscript := getCanonicalTranscript(variantAnnotation)
	// Taken from GN response (default)
	event.Chromosome = resolveChromosome(variantAnnotation, genomicLocation) // annotationUtil.resolveChromosome(gnResponse, mRecord)
	event.StartPosition = resolveStart(variantAnnotation, genomicLocation)
	event.EndPosition = resolveEnd(variantAnnotation, genomicLocation) // annotationUtil.resolveEnd(gnResponse, mRecord)
	event.Strand = resolveStrandSign(variantAnnotation)                // annotationUtil.resolveStrandSign(gnResponse, mRecord)
	event.NcbiBuild = resolveAssemblyName(variantAnnotation)           // annotationUtil.resolveAssemblyName(gnResponse, mRecord)
	event.HugoSymbol = resolveHugoSymbol(canonicalTranscript)
	event.EntrezGeneId = resolveEntrezGeneId(canonicalTranscript)
	event.VariantClassification = *canonicalTranscript.VariantClassification // annotationUtil.resolveVariantClassification(gnResponse, canonicalTranscript, mRecord)
	event.VariantType = resolveVariantType(variantAnnotation)                // annotationUtil.resolveVariantType(gnResponse)
	event.DbsnpRs = resolveDbSnpRs(variantAnnotation)                        // 	annotationUtil.resolveDbSnpRs(gnResponse, mRecord)
	event.Hgvsc = resolveHgvsc(canonicalTranscript)                          // annotationUtil.resolveHgvsc(canonicalTranscript)
	event.Hgvsp = resolveHgvsp(canonicalTranscript)                          // annotationUtil.resolveHgvsp(canonicalTranscript)
	event.HgvspShort = resolveHgvspShort(canonicalTranscript)                // annotationUtil.resolveHgvspShort(canonicalTranscript)
	event.TranscriptId = resolveTranscriptId(canonicalTranscript)            // annotationUtil.resolveTranscriptId(canonicalTranscript)
	event.Refseq = resolveRefSeq(canonicalTranscript)                        // annotationUtil.resolveRefSeq(canonicalTranscript)
	event.Codons = resolveCodonChange(canonicalTranscript)
	event.Consequence = resolveConsequence(canonicalTranscript)
	event.ProteinPosition = resolveProteinPosition(canonicalTranscript)
	event.ExonNumber = resolveExon(canonicalTranscript)

	// Taken from GN response (Polyphen)
	event.PolyphenPrediction = resolvePolyphenPrediction(canonicalTranscript) // annotationUtil.resolvePolyphenPrediction(canonicalTranscript)
	event.PolyphenScore = resolvePolyphenScore(canonicalTranscript)           // annotationUtil.resolvePolyphenScore(canonicalTranscript)

	// Taken from GN response (SIFT)
	event.SiftPrediction = resolveSiftPrediction(canonicalTranscript) // annotationUtil.resolveSiftPrediction(canonicalTranscript)
	event.SiftScore = resolveSiftScore(canonicalTranscript)           // annotationUtil.resolveSiftScore(canonicalTranscript)

	// ======================================
	// separate special function needed
	event.ReferenceAllele, event.TumorSeqAllele1, event.TumorSeqAllele2 = resolveRefAndTumorSeqAlleles(variantAnnotation, *event, "true")
	// ======================================
	event.AnnotationStatus = "SUCCESS"
	return nil
}
