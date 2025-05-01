package tempo_databricks_gateway

import (
    "context"
    "fmt"
    "os"
    "strconv"
    
    gnapi "github.com/genome-nexus/genome-nexus-go-api-client/genome-nexus-public-api"
    tt "github.mskcc.org/cdsi/cdsi-protobuf/tempo/generated/v1/go"
)

// TODO: add back in hardcoded values and/or values TEMPO must provide (currently blank)
func CreateAnnotatedEvent(variantAnnotation gnapi.VariantAnnotation, genomicLocation gnapi.GenomicLocation, event tt.Event) *tt.Event {
	if ! *variantAnnotation.SuccessfullyAnnotated {
		event.AnnotationStatus = "FAILED"
		return &event
	}
	canonicalTranscript := GetCanonicalTranscript(variantAnnotation)
	annotatedEvent := &tt.Event{}
	// Taken from GN response (default)
	annotatedEvent.Chromosome = ResolveChromosome(variantAnnotation, genomicLocation) // annotationUtil.resolveChromosome(gnResponse, mRecord)
	annotatedEvent.StartPosition = ResolveStart(variantAnnotation, genomicLocation)
	annotatedEvent.EndPosition = ResolveEnd(variantAnnotation, genomicLocation) // annotationUtil.resolveEnd(gnResponse, mRecord)
	annotatedEvent.Strand =	ResolveStrandSign(variantAnnotation) // annotationUtil.resolveStrandSign(gnResponse, mRecord)
	annotatedEvent.NcbiBuild = ResolveAssemblyName(variantAnnotation) // annotationUtil.resolveAssemblyName(gnResponse, mRecord)
	annotatedEvent.HugoSymbol = ResolveHugoSymbol(canonicalTranscript)
	annotatedEvent.EntrezGeneId = ResolveEntrezGeneId(canonicalTranscript)
	annotatedEvent.VariantClassification = *canonicalTranscript.VariantClassification // annotationUtil.resolveVariantClassification(gnResponse, canonicalTranscript, mRecord)
	annotatedEvent.VariantType = ResolveVariantType(variantAnnotation) // annotationUtil.resolveVariantType(gnResponse)
	annotatedEvent.DbsnpRs = ResolveDbSnpRs(variantAnnotation) // 	annotationUtil.resolveDbSnpRs(gnResponse, mRecord)
	annotatedEvent.Hgvsc = ResolveHgvsc(canonicalTranscript) // annotationUtil.resolveHgvsc(canonicalTranscript)
	annotatedEvent.Hgvsp = ResolveHgvsp(canonicalTranscript) // annotationUtil.resolveHgvsp(canonicalTranscript)
	annotatedEvent.HgvspShort = ResolveHgvspShort(canonicalTranscript) // annotationUtil.resolveHgvspShort(canonicalTranscript)
	annotatedEvent.TranscriptId = ResolveTranscriptId(canonicalTranscript) // annotationUtil.resolveTranscriptId(canonicalTranscript)
	annotatedEvent.Refseq = ResolveRefSeq(canonicalTranscript) // annotationUtil.resolveRefSeq(canonicalTranscript)
	annotatedEvent.Codons = ResolveCodonChange(canonicalTranscript)
	annotatedEvent.Consequence = ResolveConsequence(canonicalTranscript)
	annotatedEvent.ProteinPosition = ResolveProteinPosition(canonicalTranscript)
	annotatedEvent.ExonNumber = ResolveExon(canonicalTranscript)

	// Taken from GN response (Polyphen)
	annotatedEvent.PolyphenPrediction = ResolvePolyphenPrediction(canonicalTranscript) // annotationUtil.resolvePolyphenPrediction(canonicalTranscript)
	annotatedEvent.PolyphenScore = ResolvePolyphenScore(canonicalTranscript) // annotationUtil.resolvePolyphenScore(canonicalTranscript)

	// Taken from GN response (SIFT)
	annotatedEvent.SiftPrediction = ResolveSiftPrediction(canonicalTranscript) // annotationUtil.resolveSiftPrediction(canonicalTranscript)
	annotatedEvent.SiftScore = ResolveSiftScore(canonicalTranscript) // annotationUtil.resolveSiftScore(canonicalTranscript)

  	// ======================================
  	// separate special function needed
	annotatedEvent.ReferenceAllele, annotatedEvent.TumorSeqAllele1, annotatedEvent.TumorSeqAllele2 = ResolveRefAndTumorSeqAlleles(variantAnnotation, event, "true")
	// ======================================
	annotatedEvent.AnnotationStatus = "SUCCESS"
	return annotatedEvent
}

func AnnotateTempoMessageEvents(tm tt.TempoMessage) []*tt.Event {
        // prepare genomic locations for annotation
	genomicLocations := make([]gnapi.GenomicLocation, 0)
        for _, a := range tm.Events {
		loc := GetGenomicLocation(a)
          	genomicLocations = append(genomicLocations, loc)
        }

	// send request to genome nexus
        //fmt.Printf("About to annotate sample: %s with %s events", tm.CmoSampleId, len(genomicLocations))
        variantAnnotations := GetVariantAnnotations(genomicLocations)
        if len(variantAnnotations) != len(genomicLocations) {
        	fmt.Println("Returned number of variant annotations does not match number of genomic locations submitted")
        }

        // construct a list of annotated events using genome nexus response
        fmt.Printf("Now constructing full TempoMessage/Annotated Record for putting into S3 for sample: %s \n" , tm.CmoSampleId)
        annotatedEvents := []*tt.Event{}
        for i, variantAnnotation := range variantAnnotations {
          annotatedEvents = append(annotatedEvents, CreateAnnotatedEvent(variantAnnotation, genomicLocations[i], *tm.Events[i]))
        }
        return annotatedEvents
}

func GetVariantAnnotations(genomicLocations []gnapi.GenomicLocation) []gnapi.VariantAnnotation {
	var token = ""
  	fields := make([]string, 0)
	fields = append(fields, "annotation_summary")
	// TODO: confirm this should be mskcc and not uniprot
	var isoformOverrideSource = "mskcc"
 	configuration := gnapi.NewConfiguration()
 	apiClient := gnapi.NewAPIClient(configuration)

	ctx := context.WithValue(context.Background(), gnapi.ContextServerIndex, 1)
  	x := apiClient.AnnotationControllerApi.FetchVariantAnnotationByGenomicLocationPOST(ctx).GenomicLocations(genomicLocations)
	x = x.Fields(fields)
	x = x.IsoformOverrideSource(isoformOverrideSource)
	x = x.Token(token)

 	variantAnnotations, r, err := x.Execute()
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error when calling `AnnotationControllerApi.FetchVariantAnnotationByGenomicLocationGET``: %v\n", err)
		fmt.Fprintf(os.Stderr, "Full HTTP response: %v\n", r)
	}
	for i := 0; i < len(variantAnnotations); i++ {
		if !*variantAnnotations[i].SuccessfullyAnnotated {
			fmt.Println("ERROR: Failed to annotate: ", variantAnnotations[i].Variant)
		}
	}
	return variantAnnotations
}

func GetGenomicLocation(e *tt.Event) gnapi.GenomicLocation {
	start,_ := strconv.ParseInt(e.StartPosition,10,32)
	end,_ := strconv.ParseInt(e.StartPosition,10,32)
  
 	gloc := gnapi.NewGenomicLocation(e.Chromosome,
                                      int32(start),
                                      int32(end),
                                      e.ReferenceAllele,
                                      ResolveTumorSeqAlleleFromInput(e.ReferenceAllele, e.TumorSeqAllele1, e.TumorSeqAllele2))
	// Uncomment to get the genomic location being sent to Genome Nexus
	// fmt.Println(*gloc)
  	return *gloc
}
