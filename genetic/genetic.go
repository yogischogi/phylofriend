// Package genetic contains functions for calculating genetic distances from Y-STR values.
package genetic

import (
	"bytes"
	"errors"
	"fmt"
	"math"
	"strconv"
)

const (
	// MaxMarkers denotes the number of Y-STR markers available to the program.
	// The array that contains the markers holds additional values for
	// DYS464. So the actual array size is MaxMarkers + NDYS464ext.
	MaxMarkers = 496

	// NDYS464ext denotes the number of extra values for the DYS464 marker.
	// 98.5% of all people do not have more than four values at the DYS464 marker
	// (http://www.isogg.org/wiki/DYS_464).
	// Currently many people use spreadsheets which support only four
	// values for DYS464. To stay compatible Phylofriend also uses
	// 4 markers at the standard DYS464 position and adds the extra
	// markers at the end.
	NDYS464ext = 4

	// Marker positions:

	// DYS389i is a marker that is included in DYS389ii.
	DYS389i  = 9
	DYS389ii = 11

	// The program uses the infinite allele mutation model
	// for these palindromic markers.
	DYS464start    = 21
	DYS464end      = 24
	DYS464extStart = MaxMarkers
	DYS464extEnd   = MaxMarkers + NDYS464ext - 1
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
	// Label is used as this person's label.
	// A label must be exactly 10 characters long and may
	// only contain 8-bit characters. This is to be
	// compatible with PHYLIP and the Newick tree format.
	Label    string
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
		Label:       "__________",
		Ancestor:    "",
		Origin:      "",
		YstrMarkers: p.YstrMarkers}
}

// markerSet returns a person who's marker set has been reduced
// to the specified number of values.
// If the person has not been tested for all markers isComplete = false.
func (p *Person) markerSet(nMarkers int) (person *Person, isComplete bool) {
	person = new(Person)
	*person = *p
	isComplete = true
	// Delete all marker values outside the specified set.
	for i := nMarkers; i < MaxMarkers; i++ {
		person.YstrMarkers[i] = 0
	}
	// Check if the marker set is complete
	for i := 0; i < nMarkers; i++ {
		if person.YstrMarkers[i] == 0 &&
			// The palindromic marker DYS464 is not counted here.
			// In all cases I know, at least four values are reported
			// for DYS464, but in theory there could be less. So I
			// ignore this marker here, just in case.
			(i <= DYS464start || i >= DYS464end) {
			isComplete = false
			break
		}
	}
	return person, isComplete
}

// YstrMarkers contains the values for the Y-STR markers.
// The first 111 markers are in Family Tree DNA order.
// The detailed layout is defined by YstrMarkerTable.
type YstrMarkers [MaxMarkers + NDYS464ext]float64

func (y *YstrMarkers) String() string {
	var buffer bytes.Buffer
	for i, _ := range YstrMarkerTable {
		text := fmt.Sprintf("%s: %.f, ", YstrMarkerTable[i].InternalName, y[i])
		buffer.WriteString(text)
	}
	buffer.WriteString("\n")
	return buffer.String()
}

// DefaultMutationRates() returns a structure of YstrMarkers,
// where all values are set to 1.
// This makes mutation counting the default behaviour.
func DefaultMutationRates() YstrMarkers {
	var result YstrMarkers
	for i := 0; i < len(result); i++ {
		result[i] = 1
	}
	return result
}

// distancesMarkerCount counts the difference between every single marker provided.
// It returns a slice of distances containing the difference for every marker.
// nCompared is the number of markers compared. If one of the markers is 0, it
// is not counted. Thus people who have tested for different marker sets can
// be compared. This method is very simple and applies no special methods for
// palindromic markers.
func distancesMarkerCount(ystr1, ystr2 YstrMarkers) (distances []float64, nCompared int) {
	nCompared = 0
	distances = make([]float64, MaxMarkers)
	for i := 0; i < MaxMarkers; i++ {
		if ystr1[i] != 0 && ystr2[i] != 0 {
			distances[i] = math.Abs(ystr1[i] - ystr2[i])
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
	// nCompared is the number of markers that are actually compared.
	// We compare only those marker for which the results of two persons
	// and the mutation rate exist.
	var nCompared = 0

	// stepwise calculates the genetic distance for one marker of
	// two persons using the stepwise mutation model
	// (http://nitro.biosci.arizona.edu/ftDNA/models.html).
	var stepwise = func(marker1, marker2, mutationRate float64) (distance float64) {
		if marker1 != 0 && marker2 != 0 && mutationRate != 0 {
			distance = math.Abs(marker1-marker2) / mutationRate
			nCompared++
		}
		return distance
	}

	// infinite calculates the genetic distance for one marker of
	// two persons using the infinite allelles mutation model
	// (http://nitro.biosci.arizona.edu/ftDNA/models.html).
	var infinite = func(marker1, marker2, mutationRate float64) (distance float64) {
		if marker1 != 0 && marker2 != 0 && mutationRate != 0 {
			if marker1 != marker2 {
				distance = 1 / mutationRate
			} else {
				distance = 0
			}
			nCompared++
		}
		return distance
	}

	// Check if values for special markers exist.
	DYS389exists := false
	if ystr1[DYS389i] != 0 &&
		ystr2[DYS389i] != 0 &&
		ystr1[DYS389ii] != 0 &&
		ystr2[DYS389ii] != 0 {
		DYS389exists = true
	}
	YCAIIexists := false
	if ystr1[YCAIIa] != 0 &&
		ystr2[YCAIIa] != 0 &&
		ystr1[YCAIIb] != 0 &&
		ystr2[YCAIIb] != 0 {
		YCAIIexists = true
	}

	// Calculate distance for every single marker.
	distances := make([]float64, MaxMarkers)
	for i := 0; i < DYS389ii; i++ {
		distances[i] = stepwise(ystr1[i], ystr2[i], mutationRates[i])
	}
	if DYS389exists && mutationRates[DYS389ii] != 0 {
		distances[DYS389ii] = distanceDYS389ii(ystr1[DYS389i], ystr1[DYS389ii], ystr2[DYS389i], ystr2[DYS389ii]) / mutationRates[DYS389ii]
		nCompared++
	}
	for i := DYS389ii + 1; i < DYS464start; i++ {
		distances[i] = stepwise(ystr1[i], ystr2[i], mutationRates[i])
	}
	// DYS464: For compatibilty reasons DYS464 is stored at different range positions.
	// So we need to put all values back together.
	values1 := concat(ystr1[DYS464start:DYS464end+1], ystr1[DYS464extStart:DYS464extEnd+1])
	values2 := concat(ystr2[DYS464start:DYS464end+1], ystr2[DYS464extStart:DYS464extEnd+1])
	if isValidPalindromic(values1, values2, mutationRates[DYS464end]) {
		distances[DYS464end] = distancePalindromic(values1, values2, mutationRates[DYS464end])
		nCompared += 4
	}
	for i := DYS464end + 1; i < YCAIIa; i++ {
		distances[i] = stepwise(ystr1[i], ystr2[i], mutationRates[i])
	}
	if YCAIIexists {
		// TODO: Better use palindromic?
		// See: http://www.dna-fingerprint.com/static/PalindromicPres.pdf
		// Currently YCAIIexists is redundant.
		distances[YCAIIa] = infinite(ystr1[YCAIIa], ystr2[YCAIIa], mutationRates[YCAIIa])
		distances[YCAIIb] = infinite(ystr1[YCAIIb], ystr2[YCAIIb], mutationRates[YCAIIb])
	}
	for i := YCAIIb + 1; i < CDYstart; i++ {
		distances[i] = stepwise(ystr1[i], ystr2[i], mutationRates[i])
	}
	if isValidPalindromic(ystr1[CDYstart:CDYend+1], ystr2[CDYstart:CDYend+1], mutationRates[CDYend]) {
		distances[CDYend] = distancePalindromic(ystr1[CDYstart:CDYend+1], ystr2[CDYstart:CDYend+1], mutationRates[CDYend])
		nCompared += 2
	}
	for i := CDYend + 1; i < DYF395S1start; i++ {
		distances[i] = stepwise(ystr1[i], ystr2[i], mutationRates[i])
	}
	if isValidPalindromic(ystr1[DYF395S1start:DYF395S1end+1], ystr2[DYF395S1start:DYF395S1end+1], mutationRates[DYF395S1end]) {
		distances[DYF395S1end] = distancePalindromic(ystr1[DYF395S1start:DYF395S1end+1], ystr2[DYF395S1start:DYF395S1end+1], mutationRates[DYF395S1end])
		nCompared += 2
	}
	for i := DYF395S1end + 1; i < DYS413start; i++ {
		distances[i] = stepwise(ystr1[i], ystr2[i], mutationRates[i])
	}
	if isValidPalindromic(ystr1[DYS413start:DYS413end+1], ystr2[DYS413start:DYS413end+1], mutationRates[DYS413end]) {
		distances[DYS413end] = distancePalindromic(ystr1[DYS413start:DYS413end+1], ystr2[DYS413start:DYS413end+1], mutationRates[DYS413end])
		nCompared += 2
	}
	for i := DYS413end + 1; i < MaxMarkers; i++ {
		distances[i] = stepwise(ystr1[i], ystr2[i], mutationRates[i])
	}
	return average(distances, nCompared)
}

// distanceDYS389ii calculates the genetic distance for the DYS389ii
// marker. This marker is a special case because it includes DYS389i
//
// The input parameters are the DYS389 values for persons a and b.
func distanceDYS389ii(aDYS389i, aDYS389ii, bDYS389i, bDYS389ii float64) float64 {
	return math.Abs((aDYS389ii - aDYS389i) - (bDYS389ii - bDYS389i))
}

// distancePalincromic calculates the genetic distance for palindromic markers
// like DYS464, CDY, DYF395S1 and DYS413.
// This calculation tries to use the same approach as FamilyTreeDNA
// (https://www.familytreedna.com/learn/y-dna-testing/y-str/infinite-allele-palindromic-markers/)
// It uses the infinite allele model and differs from
// the calculation described at genebase
// (http://www.genebase.com/learning/article/46).
//
// ystr1 and ystr2 contain the palindromic values for each person.
func distancePalindromic(ystr1, ystr2 []float64, mutationRate float64) float64 {
	if mutationRate == 0 {
		return 0
	}
	var distance float64 = 0

	// Create two lists that contain only values > 0.
	list1 := make([]float64, 0, len(ystr1))
	for _, value := range ystr1 {
		if value > 0 {
			list1 = append(list1, value)
		}
	}
	list2 := make([]float64, 0, len(ystr2))
	for _, value := range ystr2 {
		if value > 0 {
			list2 = append(list2, value)
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
	// Compare each value in set1 against values of set2.
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
	return distance / mutationRate
}

// isValidPalindromic tests if the palindromic marker specified
// for two samples ystr1 and ystr2 is valid for further processing.
// The marker is considered valid if both ystr1 and ystr2 contain
// at least one value > 0 and the mutation rate is also > 0.
func isValidPalindromic(ystr1, ystr2 []float64, mutationRate float64) bool {
	if mutationRate == 0 {
		return false
	}
	isValid1 := false
	for _, value := range ystr1 {
		if value > 0 {
			isValid1 = true
			break
		}
	}
	isValid2 := false
	for _, value := range ystr2 {
		if value > 0 {
			isValid2 = true
			break
		}
	}
	return isValid1 && isValid2
}

// concat concatenates to slices and returns the result
// without any changes to the underlying arrays of s1 and s2.
func concat(s1, s2 []float64) []float64 {
	result := make([]float64, len(s1)+len(s2))
	copy(result, s1)
	copy(result[len(s1):], s2)
	return result
}

// ModalHaplotype calculates the modal haplotype for a group of persons.
// The modal value for a marker is the value with the highest occurence.
// If two values have the same frequency the lower one is chosen.
func ModalHaplotype(persons []*Person) *Person {
	modal := Person{
		ID:       "modal",
		Name:     "modal",
		Label:    "_____modal",
		Ancestor: "modal",
		Origin:   "modal",
	}
	// Calculate modal value for each marker.
	for marker := 0; marker < MaxMarkers; marker++ {
		// cMarkers maps marker values to the count of that value.
		cMarkers := make(map[float64]int)
		// Count marker values.
		for _, person := range persons {
			markerValue := person.YstrMarkers[marker]
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
		modal.YstrMarkers[marker] = modalValue
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

// Anonymize anonymizes the persons data and replaces the label
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
		result[i].Label = string(underscores) + nameNo
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
func ReduceToMarkerSet(persons []*Person, nMarkers int) ([]*Person, error) {
	result := make([]*Person, 0, len(persons))
	for _, p := range persons {
		next, isComplete := p.markerSet(nMarkers)
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
