// Package genetic contains functions for calculating genetic distances from Y-STR values.
package genetic

import (
	"errors"
	"math"
	"strconv"
)

const (
	// Nmarkers denotes the number of Y-STR markers available to the program.
	Nmarkers = 111

	// NDYS464ext denotes the number of extra values for the DYS464 marker.
	// 98.5% of all people do not have more than four values at the DYS464 marker
	// (http://www.isogg.org/wiki/DYS_464).
	// Currently many people use spreadsheets which support only four
	// values for DYS464. To stay compatible Phylofriend also uses
	// 4 markers at the standard DYS464 position and adds the extra
	// markers at the end.
	NDYS464ext = 3

	// Marker positions:

	// DYS389i is a marker that is included in DYS389ii.
	DYS389i  = 9
	DYS389ii = 11

	// The program uses the infinite allele mutation model
	// for these palindromic markers.
	DYS464start    = 21
	DYS464end      = 24
	DYS464extStart = Nmarkers
	DYS464extEnd   = Nmarkers + NDYS464ext - 1
	CDYstart       = 33
	CDYend         = 34
	DYF395S1start  = 39
	DYF395S1end    = 40
	DYS413start    = 48
	DYS413end      = 49

	// For YCAII the infinite allele mutation model is used.
	YCAIIa = 27
	YCAIIb = 28
)

// DistanceFunc is a function to calculate the genetic distance
// between two sets of Y-STR markers.
// The first two parameters are the Y-STR markers of the persons
// to compare. The last parameter contains the mutation rates
// used for comparison.
type DistanceFunc func(YstrMarkers, YstrMarkers, YstrMarkers) float64

// Person resembles a person with a set of Y-STR markers.
type Person struct {
	ID   string
	Name string
	// Name10 is a name that is exactly 10 characters long
	// and contains no Unicode characters. This is used to
	// be compatible with PHYLIP and the Newick tree format.
	Name10   string
	Ancestor string
	Origin   string
	YstrMarkers
}

// anonymize deletes personal data with the exception of the
// Y-STR values.
func (p *Person) anonymize() *Person {
	return &Person{
		ID:          "",
		Name:        "",
		Name10:      "__________",
		Ancestor:    "",
		Origin:      "",
		YstrMarkers: p.YstrMarkers}
}

// markerSet returns a person who's marker set has been reduced
// to the specified number of values.
// If the person has not been tested for all markers isComplete = false.
func (p *Person) markerSet(nValues int) (person *Person, isComplete bool) {
	person = new(Person)
	*person = *p
	isComplete = true
	// Delete all marker values outside the specified set.
	for i := nValues; i < Nmarkers; i++ {
		person.YstrMarkers.SetValue(i, 0)
	}
	// Check if the marker set is complete
	for i := 0; i < nValues; i++ {
		if person.YstrMarkers.Value(i) == 0 &&
			(i <= DYS464start || i >= DYS464end) {
			isComplete = false
			break
		}
	}
	return person, isComplete
}

// YstrMarkers contains values for 111 Y-STR markers in
// FamilyTreeDNA order (2014-06-25).
// For more information on markers take a look at
// http://en.wikipedia.org/wiki/List_of_Y-STR_markers
// or http://www.isogg.org/markersig.htm.
type YstrMarkers struct {
	DYS393,
	DYS390,
	DYS19,
	DYS391,
	DYS385a,
	DYS385b,
	DYS426,
	DYS388,
	DYS439,
	DYS389i,
	DYS392,
	DYS389ii,
	DYS458,
	DYS459a,
	DYS459b,
	DYS455,
	DYS454,
	DYS447,
	DYS437,
	DYS448,
	DYS449,
	DYS464a,
	DYS464b,
	DYS464c,
	DYS464d,
	DYS460,
	Y_GATA_H4,
	YCAIIa,
	YCAIIb,
	DYS456,
	DYS607,
	DYS576,
	DYS570,
	CDYa,
	CDYb,
	DYS442,
	DYS438,
	DYS531,
	DYS578,
	DYF395S1a,
	DYF395S1b,
	DYS590,
	DYS537,
	DYS641,
	DYS472,
	DYF406S1,
	DYS511,
	DYS425,
	DYS413a,
	DYS413b,
	DYS557,
	DYS594,
	DYS436,
	DYS490,
	DYS534,
	DYS450,
	DYS444,
	DYS481,
	DYS520,
	DYS446,
	DYS617,
	DYS568,
	DYS487,
	DYS572,
	DYS640,
	DYS492,
	DYS565,
	DYS710,
	DYS485,
	DYS632,
	DYS495,
	DYS540,
	DYS714,
	DYS716,
	DYS717,
	DYS505,
	DYS556,
	DYS549,
	DYS589,
	DYS522,
	DYS494,
	DYS533,
	DYS636,
	DYS575,
	DYS638,
	DYS462,
	DYS452,
	DYS445,
	Y_GATA_A10,
	DYS463,
	DYS441,
	Y_GGAAT_1B07,
	DYS525,
	DYS712,
	DYS593,
	DYS650,
	DYS532,
	DYS715,
	DYS504,
	DYS513,
	DYS561,
	DYS552,
	DYS726,
	DYS635,
	DYS587,
	DYS643,
	DYS497,
	DYS510,
	DYS434,
	DYS461,
	DYS435,
	DYS464e,
	DYS464f,
	DYS464g float64
}

// addressOf is a helper function for Value and SetValue to
// allow array like access to the YstrMarkers struct.
func (m *YstrMarkers) addressOf(i int) *float64 {
	markers := [Nmarkers + NDYS464ext]*float64{
		&m.DYS393,
		&m.DYS390,
		&m.DYS19,
		&m.DYS391,
		&m.DYS385a,
		&m.DYS385b,
		&m.DYS426,
		&m.DYS388,
		&m.DYS439,
		&m.DYS389i,
		&m.DYS392,
		&m.DYS389ii,
		&m.DYS458,
		&m.DYS459a,
		&m.DYS459b,
		&m.DYS455,
		&m.DYS454,
		&m.DYS447,
		&m.DYS437,
		&m.DYS448,
		&m.DYS449,
		&m.DYS464a,
		&m.DYS464b,
		&m.DYS464c,
		&m.DYS464d,
		&m.DYS460,
		&m.Y_GATA_H4,
		&m.YCAIIa,
		&m.YCAIIb,
		&m.DYS456,
		&m.DYS607,
		&m.DYS576,
		&m.DYS570,
		&m.CDYa,
		&m.CDYb,
		&m.DYS442,
		&m.DYS438,
		&m.DYS531,
		&m.DYS578,
		&m.DYF395S1a,
		&m.DYF395S1b,
		&m.DYS590,
		&m.DYS537,
		&m.DYS641,
		&m.DYS472,
		&m.DYF406S1,
		&m.DYS511,
		&m.DYS425,
		&m.DYS413a,
		&m.DYS413b,
		&m.DYS557,
		&m.DYS594,
		&m.DYS436,
		&m.DYS490,
		&m.DYS534,
		&m.DYS450,
		&m.DYS444,
		&m.DYS481,
		&m.DYS520,
		&m.DYS446,
		&m.DYS617,
		&m.DYS568,
		&m.DYS487,
		&m.DYS572,
		&m.DYS640,
		&m.DYS492,
		&m.DYS565,
		&m.DYS710,
		&m.DYS485,
		&m.DYS632,
		&m.DYS495,
		&m.DYS540,
		&m.DYS714,
		&m.DYS716,
		&m.DYS717,
		&m.DYS505,
		&m.DYS556,
		&m.DYS549,
		&m.DYS589,
		&m.DYS522,
		&m.DYS494,
		&m.DYS533,
		&m.DYS636,
		&m.DYS575,
		&m.DYS638,
		&m.DYS462,
		&m.DYS452,
		&m.DYS445,
		&m.Y_GATA_A10,
		&m.DYS463,
		&m.DYS441,
		&m.Y_GGAAT_1B07,
		&m.DYS525,
		&m.DYS712,
		&m.DYS593,
		&m.DYS650,
		&m.DYS532,
		&m.DYS715,
		&m.DYS504,
		&m.DYS513,
		&m.DYS561,
		&m.DYS552,
		&m.DYS726,
		&m.DYS635,
		&m.DYS587,
		&m.DYS643,
		&m.DYS497,
		&m.DYS510,
		&m.DYS434,
		&m.DYS461,
		&m.DYS435,
		&m.DYS464e,
		&m.DYS464f,
		&m.DYS464g}
	return markers[i]
}

// Value returns the i-th value of a set of Y-STR markers.
func (m *YstrMarkers) Value(i int) float64 {
	return *m.addressOf(i)
}

// SetValue sets the i-th value of a set of Y-STR markers.
func (m *YstrMarkers) SetValue(i int, value float64) {
	*m.addressOf(i) = value
}

// MutationRates is a set of mutation rates that are
// used as default values.
// Here 37 marker average values are used. They can be found in
// http://www.scirp.org/journal/PaperInformation.aspx?paperID=19567#.U-IUh9a3KR8.
var MutationRates = YstrMarkers{
	DYS393:    0.00243,
	DYS390:    0.00243,
	DYS19:     0.00243,
	DYS391:    0.00243,
	DYS385a:   0.00243,
	DYS385b:   0.00243,
	DYS426:    0.00243,
	DYS388:    0.00243,
	DYS439:    0.00243,
	DYS389i:   0.00243,
	DYS392:    0.00243,
	DYS389ii:  0.00243,
	DYS458:    0.00243,
	DYS459a:   0.00243,
	DYS459b:   0.00243,
	DYS455:    0.00243,
	DYS454:    0.00243,
	DYS447:    0.00243,
	DYS437:    0.00243,
	DYS448:    0.00243,
	DYS449:    0.00243,
	DYS464a:   0.00243,
	DYS464b:   0.00243,
	DYS464c:   0.00243,
	DYS464d:   0.00243,
	DYS460:    0.00243,
	Y_GATA_H4: 0.00243,
	YCAIIa:    0.00243,
	YCAIIb:    0.00243,
	DYS456:    0.00243,
	DYS607:    0.00243,
	DYS576:    0.00243,
	DYS570:    0.00243,
	CDYa:      0.00243,
	CDYb:      0.00243,
	DYS442:    0.00243,
	DYS438:    0.00243}

// distancesMarkerCount counts the difference between every single marker provided.
// It returns a slice of distances containing the difference for every marker.
// nCompared is the number of markers compared. If one of the markers is 0, it
// is not counted. Thus people who have tested for different marker sets can
// be compared. This method is very simple and applies no special methods for
// palindromic markers.
func distancesMarkerCount(ystr1, ystr2 YstrMarkers) (distances []float64, nCompared int) {
	nCompared = 0
	distances = make([]float64, Nmarkers)
	for i := 0; i < Nmarkers; i++ {
		if ystr1.Value(i) != 0 && ystr2.Value(i) != 0 {
			distances[i] = math.Abs(ystr1.Value(i) - ystr2.Value(i))
			nCompared++
		}
	}
	return distances, nCompared
}

// sum returns the sum of all values in the distances slice.
func sum(distances []float64) float64 {
	var sum float64 = 0
	for _, value := range distances {
		sum += value
	}
	return sum
}

// average calculates the average value for the distances slice.
func average(distances []float64, nCompared int) float64 {
	return sum(distances) / float64(nCompared)
}

// distancesSimpleCount returns an average genetic distance just by
// counting the marker differences.
// This should only be used for debugging purposes because the result
// is no accurate genetic distance.
// The parameter mutationRates is ignored and may be nil.
func distanceSimpleCount(ystr1, ystr2, mutationRates YstrMarkers) float64 {
	distances, nCompared := distancesMarkerCount(ystr1, ystr2)
	return average(distances, nCompared)
}

// Distance calculates the genetic distance between two sets of
// Y-STR markers.
// It uses a hybrid mutation model. For most markers the stepwise
// mutation model is used but for palindromic markers the infinite
// allele model is used. More information about mutation models
// can be found at http://nitro.biosci.arizona.edu/ftDNA/models.html.
// If one value or the mutation rate for a specific marker is
// set to 0 it is excluded from the calculation.
//
// This method may change in future versions.
func Distance(ystr1, ystr2, mutationRates YstrMarkers) float64 {
	// Check if values for special markers exist.
	DYS389exists := false
	DYS464exists := false
	YCAIIexists := false
	CDYexists := false
	DYF395S1exists := false
	DYS413exists := false

	if ystr1.Value(DYS389i) != 0 &&
		ystr2.Value(DYS389i) != 0 &&
		ystr1.Value(DYS389ii) != 0 &&
		ystr2.Value(DYS389ii) != 0 {
		DYS389exists = true
	}

	// Make sure that at least one value is reported for DYS464.
	DYS464firstExists := false
	DYS464secondExists := false
	for i := DYS464start; i <= DYS464end; i++ {
		if ystr1.Value(i) != 0 {
			DYS464firstExists = true
		}
		if ystr2.Value(i) != 0 {
			DYS464secondExists = true
		}
		if DYS464firstExists && DYS464secondExists {
			DYS464exists = true
			break
		}
	}

	if ystr1.Value(DYS464start) != 0 &&
		ystr2.Value(DYS464start) != 0 {
		DYS464exists = true
	}

	if ystr1.Value(YCAIIa) != 0 &&
		ystr2.Value(YCAIIa) != 0 &&
		ystr1.Value(YCAIIb) != 0 &&
		ystr2.Value(YCAIIb) != 0 {
		YCAIIexists = true
	}

	if ystr1.Value(CDYstart) != 0 &&
		ystr2.Value(CDYstart) != 0 {
		CDYexists = true
	}

	if ystr1.Value(DYF395S1start) != 0 &&
		ystr2.Value(DYF395S1start) != 0 {
		DYF395S1exists = true
	}

	if ystr1.Value(DYS413start) != 0 &&
		ystr2.Value(DYS413start) != 0 {
		DYS413exists = true
	}

	// Calculate distance for every single marker.
	distances := make([]float64, Nmarkers)
	// Set nCompared to 0 because we count only the values where
	// also a mutation rate exists.
	nCompared := 0
	for i := 0; i < DYS389i; i++ {
		if ystr1.Value(i) != 0 && ystr2.Value(i) != 0 && mutationRates.Value(i) != 0 {
			distances[i] = math.Abs(ystr1.Value(i)-ystr2.Value(i)) / mutationRates.Value(i)
			nCompared++
		}
	}
	if DYS389exists && mutationRates.Value(DYS389i) != 0 {
		distances[DYS389i] = distanceInfiniteAllele(ystr1, ystr2, DYS389i) / mutationRates.Value(DYS389i)
		nCompared++
	}
	for i := DYS389i + 1; i < DYS389ii; i++ {
		if ystr1.Value(i) != 0 && ystr2.Value(i) != 0 && mutationRates.Value(i) != 0 {
			distances[i] = math.Abs(ystr1.Value(i)-ystr2.Value(i)) / mutationRates.Value(i)
			nCompared++
		}
	}
	if DYS389exists && mutationRates.Value(DYS389ii) != 0 {
		distances[DYS389ii] = distanceDys389ii(ystr1, ystr2, DYS389i, DYS389ii) / mutationRates.Value(DYS389i)
		nCompared++
	}
	for i := DYS389ii + 1; i < DYS464start; i++ {
		if ystr1.Value(i) != 0 && ystr2.Value(i) != 0 && mutationRates.Value(i) != 0 {
			distances[i] = math.Abs(ystr1.Value(i)-ystr2.Value(i)) / mutationRates.Value(i)
			nCompared++
		}
	}
	if DYS464exists && mutationRates.Value(DYS464end) != 0 {
		distances[DYS464end] = distancePalindromic(ystr1, ystr2, DYS464start, DYS464end) / mutationRates.Value(DYS464end)
		nCompared += 4
	}
	for i := DYS464end + 1; i < YCAIIa; i++ {
		if ystr1.Value(i) != 0 && ystr2.Value(i) != 0 && mutationRates.Value(i) != 0 {
			distances[i] = math.Abs(ystr1.Value(i)-ystr2.Value(i)) / mutationRates.Value(i)
			nCompared++
		}
	}
	if YCAIIexists && mutationRates.Value(YCAIIa) != 0 {
		distances[YCAIIa] = distanceInfiniteAllele(ystr1, ystr2, YCAIIa) / mutationRates.Value(YCAIIa)
		nCompared++
	}
	if YCAIIexists && mutationRates.Value(YCAIIb) != 0 {
		distances[YCAIIb] = distanceInfiniteAllele(ystr1, ystr2, YCAIIb) / mutationRates.Value(YCAIIb)
		nCompared++
	}
	for i := YCAIIb + 1; i < CDYstart; i++ {
		if ystr1.Value(i) != 0 && ystr2.Value(i) != 0 && mutationRates.Value(i) != 0 {
			distances[i] = math.Abs(ystr1.Value(i)-ystr2.Value(i)) / mutationRates.Value(i)
			nCompared++
		}
	}
	if CDYexists && mutationRates.Value(CDYend) != 0 {
		distances[CDYend] = distancePalindromic(ystr1, ystr2, CDYstart, CDYend) / mutationRates.Value(CDYend)
		nCompared += 2
	}
	for i := CDYend + 1; i < DYF395S1start; i++ {
		if ystr1.Value(i) != 0 && ystr2.Value(i) != 0 && mutationRates.Value(i) != 0 {
			distances[i] = math.Abs(ystr1.Value(i)-ystr2.Value(i)) / mutationRates.Value(i)
			nCompared++
		}
	}
	if DYF395S1exists && mutationRates.Value(DYF395S1end) != 0 {
		distances[DYF395S1end] = distancePalindromic(ystr1, ystr2, DYF395S1start, DYF395S1end) / mutationRates.Value(DYF395S1end)
		nCompared += 2
	}
	for i := DYF395S1end + 1; i < DYS413start; i++ {
		if ystr1.Value(i) != 0 && ystr2.Value(i) != 0 && mutationRates.Value(i) != 0 {
			distances[i] = math.Abs(ystr1.Value(i)-ystr2.Value(i)) / mutationRates.Value(i)
			nCompared++
		}
	}
	if DYS413exists && mutationRates.Value(DYS413end) != 0 {
		distances[DYS413end] = distancePalindromic(ystr1, ystr2, DYS413start, DYS413end) / mutationRates.Value(DYS413end)
		nCompared += 2
	}
	for i := DYS413end + 1; i < Nmarkers; i++ {
		if ystr1.Value(i) != 0 && ystr2.Value(i) != 0 && mutationRates.Value(i) != 0 {
			distances[i] = math.Abs(ystr1.Value(i)-ystr2.Value(i)) / mutationRates.Value(i)
			nCompared++
		}
	}
	return average(distances, nCompared)
}

// distanceInfiniteAllele calculates the genetic distance between
// two sets of Y-STR markers at position a. It uses the infinite
// allele mutation model
// (http://nitro.biosci.arizona.edu/ftDNA/models.html)
// in which single marker differences are counted as one.
func distanceInfiniteAllele(ystr1, ystr2 YstrMarkers, a int) float64 {
	if ystr1.Value(a) != ystr2.Value(a) {
		return 1
	} else {
		return 0
	}
}

// distanceDys389ii calculates the genetic distance for the DYS389ii
// marker. This marker is a special case because DYS389i is included
// in DYS389ii (http://www.genebase.com/learning/article/46).
//
// a: Position of DYS389i
//
// b: Position of DYS389ii
func distanceDys389ii(ystr1, ystr2 YstrMarkers, a, b int) float64 {
	distii := math.Abs((ystr1.Value(b) - ystr1.Value(a)) - (ystr2.Value(b) - ystr2.Value(a)))
	return distii
}

// distancePalincromic calculates the genetic distance for palindromic markers
// like DYS464, CDY, DYF395S1 and DYS413.
// This calculation tries to use the same approach as FamilyTreeDNA
// (https://www.familytreedna.com/learn/y-dna-testing/y-str/infinite-allele-palindromic-markers/)
// It uses the infinite allele model and differs from
// the calculation described at genebase
// (http://www.genebase.com/learning/article/46).
//
// ystr1 and ystr2 contain a full set of Y-STR markers from two persons.
// a and b denote start and end position of the palindromic marker
// to be used for calculation.
func distancePalindromic(ystr1, ystr2 YstrMarkers, a, b int) float64 {
	var distance float64 = 0
	// Create two lists of values
	list1 := make([]float64, 0, b-a+1)
	list2 := make([]float64, 0, b-a+1)
	for i := a; i <= b; i++ {
		if ystr1.Value(i) > 0 {
			list1 = append(list1, ystr1.Value(i))
		}
		if ystr2.Value(i) > 0 {
			list2 = append(list2, ystr2.Value(i))
		}
	}
	// Different lengths count as 1.
	if len(list1) != len(list2) {
		distance = 1
		// Make sure set1 is smaller than set2
		if len(list1) > len(list2) {
			list1, list2 = list2, list1
		}
	}
	// Check every key in set1 against keys of set2.
	for i := 0; i < len(list1); i++ {
		isMatchingMarker := false
		for j := 0; j < len(list2); j++ {
			if list1[i] == list2[j] {
				isMatchingMarker = true
				// Set to 0 because it had a match.
				list2[j] = 0
				break
			}
		}
		if isMatchingMarker == false {
			distance++
		}
	}
	return distance
}

// ModalHaplotype calculates the modal haplotype for a group of persons.
// The modal value for a marker is the value with the highest occurence.
// If two values have the same frequency the lower one is chosen.
func ModalHaplotype(persons []*Person) *Person {
	modal := Person{
		ID:       "modal",
		Name:     "modal",
		Name10:   "_____modal",
		Ancestor: "modal",
		Origin:   "modal",
	}
	// Calculate modal value for each marker.
	for marker := 0; marker < Nmarkers; marker++ {
		// cMarkers maps marker values to the count of that value.
		cMarkers := make(map[float64]int)
		// Count marker values.
		for _, person := range persons {
			markerValue := person.YstrMarkers.Value(marker)
			if markerValue > 0 {
				cMarkers[markerValue] += 1
			}
		}
		// Find modal value.
		max := 0
		modalValue := 0.0
		// Determine the marker with the highest occurence.
		for value, count := range cMarkers {
			if count > max {
				modalValue = value
				max = count
			}
		}
		// If two markers have the same frequency choose the lower one.
		for value, count := range cMarkers {
			if count == max && value < modalValue {
				modalValue = value
			}
		}
		modal.YstrMarkers.SetValue(marker, modalValue)
	}
	return &modal
}

// DistanceMatrix is a matrix of genetic distances for a list of persons.
// Distance matrices are used as input for phylogenetic tree software
// like PHYLIP (https://en.wikipedia.org/wiki/PHYLIP).
type DistanceMatrix struct {
	Size   int
	Values [][]float64
}

// NewDistanceMatrix creates a genetic distance matrix for a list of persons.
func NewDistanceMatrix(
	persons []*Person,
	mutationRates YstrMarkers,
	distance DistanceFunc,
) *DistanceMatrix {
	matrix := new(DistanceMatrix)
	matrix.Size = len(persons)

	// Allocate space.
	matrix.Values = make([][]float64, matrix.Size)
	for line := 0; line < matrix.Size; line++ {
		matrix.Values[line] = make([]float64, matrix.Size)
	}

	// Calculate genetic distances for the upper right triangle.
	for i := 0; i < matrix.Size; i++ {
		for j := i; j < matrix.Size; j++ {
			matrix.Values[i][j] = distance(persons[i].YstrMarkers, persons[j].YstrMarkers, mutationRates)
		}
	}

	// Calculate genetic distances for the lower left triangle.
	for i := 1; i < matrix.Size; i++ {
		for j := 0; j < i; j++ {
			matrix.Values[i][j] = matrix.Values[j][i]
		}
	}
	return matrix
}

// Years returns a new Distance matrix that contains the distances in years units.
// The entries are just multiplied by generationDistance and calibrationFactor.
//
// generationDistance is usually something from 20 to 35. 25 and 30 seem to
// be most wildly used.
//
// calibrationFactor is used to calibrate the data to fit historical events.
func (dm *DistanceMatrix) Years(generationDistance, calibrationFactor float64) *DistanceMatrix {
	factor := generationDistance * calibrationFactor
	result := new(DistanceMatrix)
	result.Size = dm.Size

	// Allocate space
	result.Values = make([][]float64, result.Size)
	for row := 0; row < result.Size; row++ {
		result.Values[row] = make([]float64, result.Size)
	}

	// Multiply values by calibration factor.
	for i := 0; i < result.Size; i++ {
		for j := 0; j < result.Size; j++ {
			result.Values[i][j] = math.Trunc(factor * dm.Values[i][j])
		}
	}
	return result
}

// Anonymize anonymizes the persons data and replaces the name
// by a number.
func Anonymize(persons []*Person) []*Person {
	result := make([]*Person, len(persons))
	for i, p := range persons {
		result[i] = p.anonymize()
		// Create name with the exact lenght of 10
		// so that it satisfies the Newick tree specification.
		nameNo := strconv.Itoa(i + 1)
		underscores := make([]byte, 10-len(nameNo))
		for j, _ := range underscores {
			underscores[j] = '_'
		}
		result[i].Name10 = string(underscores) + nameNo
	}
	return result
}

// Reduce reduces the number of persons by dividing through factor.
// The method returns every factor-th person from the
// input data set and returns an error if the reduced data
// set contains less than two persons.
// This fuction is intended for large data sets that take a
// long time to process.
func Reduce(persons []*Person, factor int) ([]*Person, error) {
	result := make([]*Person, len(persons)/factor)
	if len(result) < 2 {
		return result, errors.New("reduction too large")
	}
	for i, _ := range result {
		result[i] = persons[i*factor]
	}
	return result, nil
}

// ReduceToMarkerSet reduces a slice of persons to only those
// persons who have tested at least for the specified number of markers.
// Additional markers are set to 0.
// If the result contains less than two persons this function returns an error.
func ReduceToMarkerSet(persons []*Person, nValues int) ([]*Person, error) {
	result := make([]*Person, 0, len(persons))
	for _, p := range persons {
		next, isComplete := p.markerSet(nValues)
		if isComplete == true {
			result = append(result, next)
		}
	}
	if len(result) < 2 {
		return result, errors.New("not enough persons who have tested for so many markers")
	}
	return result, nil
}

// Average calculates the average and the standard deviation
// for a slice of values.
// The standard deviation is calculated by using the following
// estimate: s = Sqrt(1/(N-1) * Sum((x_i - m)^2))
func Average(values []float64) (m, s float64, err error) {
	if len(values) < 2 {
		return m, s, errors.New("too few values")
	}
	N := float64(len(values))
	m = sum(values) / N
	d := 0.0
	for i, _ := range values {
		d += (values[i] - m) * (values[i] - m)
	}
	s = math.Sqrt(d / (N - 1))
	return m, s, err
}
