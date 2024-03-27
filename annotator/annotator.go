package annotator

import (
	"context"
	"fmt"
	"os"

	pubgn "github.com/mskcc/mskprotobuf/genome-nexus-public-api"
)

// TODO add summary statistics
// Java call: annotator.getAnnotatedRecordsUsingPOST(summaryStatistics, records, "mskcc", true, postIntervalSize, reannotate);

func GetVariantAnnotations(genomicLocations []pubgn.GenomicLocation) []pubgn.VariantAnnotation {
	var token = ""
	fields := make([]string, 0)
	fields = append(fields, "annotation_summary")
	var isoformOverrideSource = ""
	configuration := pubgn.NewConfiguration()
	apiClient := pubgn.NewAPIClient(configuration)
	x := apiClient.AnnotationControllerApi.FetchVariantAnnotationByGenomicLocationPOST(context.Background()).GenomicLocations(genomicLocations)
	x = x.Fields(fields)
	x = x.IsoformOverrideSource(isoformOverrideSource)
	x = x.Token(token)
	variantAnnotations, r, err := x.Execute()
	if err != nil {
		fmt.Fprintf(os.Stderr, "Error when calling `AnnotationControllerApi.FetchVariantAnnotationByGenomicLocationGET``: %v\n", err)
		fmt.Fprintf(os.Stderr, "Full HTTP response: %v\n", r)
	}
	fmt.Println("sigh \n\n\n\n")
	fmt.Println(r.Request)
	fmt.Println(r.Header)
	// https://github.com/mskcc/go-mpath-proto/blob/16d6085cee926ad85fc29c6c62319ee3697edcf2/genome-nexus-public-api/api_annotation_controller.go#L218
	// response from `FetchVariantAnnotationByGenomicLocationPOST`: []VariantAnnotation, *http.Response, error
	// VariantAnnotation defined here: https://github.com/mskcc/go-mpath-proto/blob/16d6085cee926ad85fc29c6c62319ee3697edcf2/genome-nexus-internal-api/model_variant_annotation.go#L21
	for i := 0; i < len(variantAnnotations); i++ {
		if !*variantAnnotations[i].SuccessfullyAnnotated {
			fmt.Println("ERROR: Failed to annotate: ", variantAnnotations[i].Variant)
		} else {
			fmt.Println("Success: ", *variantAnnotations[i].SuccessfullyAnnotated)
		}
	}
	return variantAnnotations
}
