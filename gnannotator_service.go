package tempo_databricks_gateway

import (
    "context"
    "fmt"
    "os"
    "strconv"
    
    gnapi "github.com/genome-nexus/genome-nexus-go-api-client/genome-nexus-public-api"
    tt "github.mskcc.org/cdsi/cdsi-protobuf/tempo/generated/v1/go"
)

func AnnotateTempoMessageEvents(tm *tt.TempoMessage) error {
        // prepare genomic locations for annotation
	genomicLocations := make([]gnapi.GenomicLocation, 0)
        for _, a := range tm.Events {
		loc := GetGenomicLocation(a)
          	genomicLocations = append(genomicLocations, loc)
        }

        variantAnnotations, err := GetVariantAnnotations(genomicLocations)
	if err != nil {
		return err
	}
        if len(variantAnnotations) != len(genomicLocations) {
		return fmt.Errorf("Returned number of variant annotations does not match genomic locations submitted")
        }

        // construct a list of annotated events using genome nexus response
        for i, variantAnnotation := range variantAnnotations {
         	mapResponseToEvent(variantAnnotation, genomicLocations[i], tm.Events[i])
        }
	return nil
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

func GetVariantAnnotations(genomicLocations []gnapi.GenomicLocation) ([]gnapi.VariantAnnotation, error) {
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
		return variantAnnotations, err
	}
	return variantAnnotations, nil
}

func mapResponseToEvent(variantAnnotation gnapi.VariantAnnotation, genomicLocation gnapi.GenomicLocation, event *tt.Event) error {
	if ! *variantAnnotation.SuccessfullyAnnotated {
		event.AnnotationStatus = "FAILED"
		return nil
	}
	canonicalTranscript := GetCanonicalTranscript(variantAnnotation)
	// Taken from GN response (default)
	event.Chromosome = ResolveChromosome(variantAnnotation, genomicLocation) // annotationUtil.resolveChromosome(gnResponse, mRecord)
	event.StartPosition = ResolveStart(variantAnnotation, genomicLocation)
	event.EndPosition = ResolveEnd(variantAnnotation, genomicLocation) // annotationUtil.resolveEnd(gnResponse, mRecord)
	event.Strand =	ResolveStrandSign(variantAnnotation) // annotationUtil.resolveStrandSign(gnResponse, mRecord)
	event.NcbiBuild = ResolveAssemblyName(variantAnnotation) // annotationUtil.resolveAssemblyName(gnResponse, mRecord)
	event.HugoSymbol = ResolveHugoSymbol(canonicalTranscript)
	event.EntrezGeneId = ResolveEntrezGeneId(canonicalTranscript)
	event.VariantClassification = *canonicalTranscript.VariantClassification // annotationUtil.resolveVariantClassification(gnResponse, canonicalTranscript, mRecord)
	event.VariantType = ResolveVariantType(variantAnnotation) // annotationUtil.resolveVariantType(gnResponse)
	event.DbsnpRs = ResolveDbSnpRs(variantAnnotation) // 	annotationUtil.resolveDbSnpRs(gnResponse, mRecord)
	event.Hgvsc = ResolveHgvsc(canonicalTranscript) // annotationUtil.resolveHgvsc(canonicalTranscript)
	event.Hgvsp = ResolveHgvsp(canonicalTranscript) // annotationUtil.resolveHgvsp(canonicalTranscript)
	event.HgvspShort = ResolveHgvspShort(canonicalTranscript) // annotationUtil.resolveHgvspShort(canonicalTranscript)
	event.TranscriptId = ResolveTranscriptId(canonicalTranscript) // annotationUtil.resolveTranscriptId(canonicalTranscript)
	event.Refseq = ResolveRefSeq(canonicalTranscript) // annotationUtil.resolveRefSeq(canonicalTranscript)
	event.Codons = ResolveCodonChange(canonicalTranscript)
	event.Consequence = ResolveConsequence(canonicalTranscript)
	event.ProteinPosition = ResolveProteinPosition(canonicalTranscript)
	event.ExonNumber = ResolveExon(canonicalTranscript)

	// Taken from GN response (Polyphen)
	event.PolyphenPrediction = ResolvePolyphenPrediction(canonicalTranscript) // annotationUtil.resolvePolyphenPrediction(canonicalTranscript)
	event.PolyphenScore = ResolvePolyphenScore(canonicalTranscript) // annotationUtil.resolvePolyphenScore(canonicalTranscript)

	// Taken from GN response (SIFT)
	event.SiftPrediction = ResolveSiftPrediction(canonicalTranscript) // annotationUtil.resolveSiftPrediction(canonicalTranscript)
	event.SiftScore = ResolveSiftScore(canonicalTranscript) // annotationUtil.resolveSiftScore(canonicalTranscript)

  	// ======================================
  	// separate special function needed
	event.ReferenceAllele, event.TumorSeqAllele1, event.TumorSeqAllele2 = ResolveRefAndTumorSeqAlleles(variantAnnotation, *event, "true")
	// ======================================
	event.AnnotationStatus = "SUCCESS"
	return nil
}

