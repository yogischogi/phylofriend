// Package genfiles provides common file operations for genetic data.
package genfiles

import (
	"bufio"
	"bytes"
	"encoding/csv"
	"encoding/json"
	"errors"
	"fmt"
	"io/ioutil"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"unicode/utf8"

	"github.com/yogischogi/phylofriend/genetic"
)

// ReadPersonsFromCSV reads persons' data from a CSV file.
// The file may contain comments or other non data rows.
// The function will try to recognize and remove them.
// The first entry of each line containing person data
// must be a unique ID.
// The Y-STR values must be in Family Tree DNA order.
// Missing values are set to 0.
//
// Some people use spreadsheets where DYS464 is stored in
// four different columns. Family Tree DNA stores it in a
// single column where the values are separated by "-".
// The function tries to recognize the storage format and
// to handle it appropriately.
//
// labelCol is the number of the colum used as a label for
// the person.
func ReadPersonsFromCSV(filename string, labelCol int) ([]*genetic.Person, error) {
	infile, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer infile.Close()

	// Read all CSV records from file.
	csvReader := csv.NewReader(infile)
	records, err := csvReader.ReadAll()
	if err != nil {
		return nil, err
	}

	// Extract lines that contain data of a sample.
	sampleRecords := make([][]string, 0, 1000)
	strIdx := 0
	for _, record := range records {
		strIdx = isSampleRecord(record)
		if strIdx > 0 {
			sampleRecords = append(sampleRecords, record)
		}
	}

	// Try to determine file format.
	// If the file format is Family Tree DNA, then DYS464
	// values are separated by a "-".
	DYS464idx := strIdx + 19
	isFTDNA := false
	for _, record := range sampleRecords {
		if strings.Contains(record[DYS464idx], "-") {
			isFTDNA = true
			break
		}
	}

	// Extract persons data from CSV records.
	persons := make([]*genetic.Person, 0, 1000)
	for _, record := range sampleRecords {
		person, err := personFromFields(record, labelCol, strIdx, isFTDNA)
		if err == nil {
			persons = append(persons, person)
		}
	}
	return persons, nil
}

// sampleRecord tries to determine if a given record is a
// valid entry of a sample containing STR data.
// If this is true, the function returns the start index
// of the STR data. Otherwise strIdx = 0.
func isSampleRecord(fields []string) (strIdx int) {
	// A valid entry must have at least one Id and 12 Y-STR values.
	if len(fields) < 13 {
		return 0
	}
	// Try to determine the start of the Y-STR values.
	// Results starts with DYS393, which is in the range 9-17.
	for i, field := range fields {
		dys393, err := strconv.ParseFloat(field, 64)
		if err == nil && dys393 >= 9 && dys393 < 17 {
			strIdx = i
			break
		}
	}
	return strIdx
}

// personFromFields creates a person from a slice of strings.
//
// The first field must contain the ID of the person.
// labelIdx is the field's index that is used for the person's Label field.
// strIdx is the index of the first STR marker value.
// isFTDNA determines if the format of the fields is Family Tree DNA like
// with palindromic values separated by "-" or not.
func personFromFields(fields []string, labelIdx, strIdx int, isFTDNA bool) (*genetic.Person, error) {
	var person genetic.Person
	var err error

	// The ID is usually the Kit number.
	person.ID = strings.TrimSpace(fields[0])
	// Return if the ID is too short.
	if len(person.ID) < 1 {
		return nil, errors.New("could not determine person ID")
	}

	person.Label = stringToLabel(strings.TrimSpace(fields[labelIdx]))
	if isFTDNA {
		person.YstrMarkers, err = extractYstrMarkersFTDNA(fields[strIdx:])
	} else {
		person.YstrMarkers, err = extractYstrMarkers(fields[strIdx:])
	}
	return &person, err
}

// extractYstrMarkers creates YstrMarkers from a slice of text fields.
// The entries must be in FamilyTreeDNA order.
func extractYstrMarkers(fields []string) (genetic.YstrMarkers, error) {
	var markers genetic.YstrMarkers

	// Trim whitespaces.
	for i, _ := range fields {
		fields[i] = strings.TrimSpace(fields[i])
	}

	// Convert string values to YstrMarkers.
	for i, strValue := range fields {
		value, err := stringToSTR(strValue)
		if err != nil {
			return markers, err
		}
		markers[i] = value
	}
	return markers, nil
}

// extractYstrMarkersFTDNA creates YstrMarkers from a slice of text fields.
// The entries must be in FamilyTreeDNA order and palindromic marker values
// must be separated by "-".
// Palindromic markers or markers with possible multiple values:
//
//   DYS19,  1 value,   idx 2,
//   DYS385, 2 values,  idx 4,
//   DYS459, 2 values, idx 12,
//   DYS464, 4 values, idx 19,
//   YCAII,  2 values, idx 22,
//   CDY,    2 values, idx 27,
//   DYF395S1, 2 values, idx 32,
//   DYS413,   2 values, idx 40.
func extractYstrMarkersFTDNA(fields []string) (genetic.YstrMarkers, error) {
	const (
		DYS385idx   = 4
		DYS459idx   = 12
		DYS464idx   = 19
		YCAIIidx    = 22
		CDYidx      = 27
		DYF395S1idx = 32
		DYS413idx   = 40
	)
	var markers genetic.YstrMarkers

	// Trim whitespaces.
	for i, _ := range fields {
		fields[i] = strings.TrimSpace(fields[i])
	}

	outIdx := 0
	for i, strValue := range fields {
		switch i {
		case DYS385idx, DYS459idx, YCAIIidx, CDYidx, DYF395S1idx, DYS413idx:
			// Extract 2 markers.
			values, err := extractMarkersFromString(strValue, 2)
			if err != nil {
				return markers, err
			}
			copy(markers[outIdx:], values)
			outIdx += 2
		case DYS464idx:
			// Extract 4 to 8 markers.
			values, err := extractMarkersFromString(strValue, 8)
			if err != nil {
				return markers, err
			}
			copy(markers[genetic.DYS464start:], values[0:4])
			copy(markers[genetic.DYS464extStart:], values[4:8])
			outIdx += 4
		default:
			// Extract one marker.
			values, err := extractMarkersFromString(strValue, 1)
			if err != nil {
				return markers, err
			}
			markers[outIdx] = values[0]
			outIdx++
		}
	}
	return markers, nil
}

// extractMarkersFromString extracts marker values from a text field
// where values are separated by "-".
// nMarkers is the maximum number of values that should be returned.
// The returned slice has always the size nMarkers. If no marker
// could be extracted the corresponding return field has the value 0.
func extractMarkersFromString(text string, nMarkers int) ([]float64, error) {
	result := make([]float64, nMarkers, nMarkers)
	fields := strings.Split(text, "-")
	count := nMarkers
	if len(fields) < count {
		count = len(fields)
	}
	for i := 0; i < count; i++ {
		value, err := stringToSTR(fields[i])
		if err != nil {
			return result, err
		}
		result[i] = value
	}
	return result, nil
}

// stringToStr converts an Y-STR value in string format to a number.
func stringToSTR(value string) (float64, error) {
	var result float64
	var err error
	if value == "" || value == "O" {
		result = 0
	} else {
		result, err = strconv.ParseFloat(value, 64)
		if err != nil {
			return result, err
		}
	}
	return result, nil
}

// ReadPersonsFromTXT reads persons' data from a text file.
// The file may contain comment lines starting with //.
// The first entry of each line is used for the person's
// Label field that is used for output in distance matrices.
// Missing Y-STR values are set to 0.
func ReadPersonsFromTXT(filename string) ([]*genetic.Person, error) {
	lines := make([]string, 0, 1000)

	// Open file
	infile, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer infile.Close()

	// Read lines
	scanner := bufio.NewScanner(infile)
	for scanner.Scan() {
		// Name must have exactly 10 characters.
		// This makes sure that empty lines are ignored
		// and comments are left out.
		if len(scanner.Text()) > 10 && scanner.Text()[0:2] != "//" {
			lines = append(lines, scanner.Text())
		}
	}
	if scanner.Err() != nil {
		return nil, scanner.Err()
	}

	// Create Persons with YstrMarkers
	persons := make([]*genetic.Person, len(lines))
	for i, _ := range lines {
		fields := strings.Fields(lines[i])
		persons[i] = new(genetic.Person)
		persons[i].Label = fields[0]
		nValues := len(fields) - 1
		if genetic.MaxMarkers < nValues {
			nValues = genetic.MaxMarkers
		}
		for j := 0; j < nValues; j++ {
			value, err := strconv.ParseFloat(fields[j+1], 64)
			if err != nil {
				return persons, err
			}
			persons[i].YstrMarkers[j] = value
		}
	}
	return persons, err
}

// ReadPersonsFromDir reads persons from the specified directory.
// All files including data must have the extension ".csv" and be
// in YFull Y-STR data format.
func ReadPersonsFromDir(dirName string) ([]*genetic.Person, error) {
	result := make([]*genetic.Person, 0, 100)
	// Get list of input files.
	infiles, err := namesWithExt(dirName, ".csv")
	if err != nil {
		return result, errors.New(fmt.Sprintf("could not read files from directory, %s\n", err))
	}
	// Read Y-STR data from input files.
	for _, infile := range infiles {
		person, err := ReadPersonFromYFull(filepath.Join(dirName, infile))
		if err != nil {
			// We do not return an error here because a
			// single invalid file should not terminate the
			// whole program.
			fmt.Printf("Error! Could not read person from file %s, %s\n", infile, err)
		} else {
			result = append(result, person)
		}
	}
	return result, nil
}

// ReadPersonFromYFull reads a person from a YFull Y-STR results file.
// The Y-Full ID is extracted from the file name and used as the persons's ID.
func ReadPersonFromYFull(filename string) (*genetic.Person, error) {
	infile, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer infile.Close()

	// Read all CSV records from file.
	csvReader := csv.NewReader(infile)
	csvReader.Comma = ';'
	records, err := csvReader.ReadAll()
	switch {
	case err != nil:
		return nil, err
	case len(records) == 0:
		return nil, errors.New(fmt.Sprintf("no data found in %s", infile))
	case len(records[0]) < 2:
		return nil, errors.New(fmt.Sprintf("invalid file format for %s", infile))
	}

	// Extract Y-STR marker values.
	var result genetic.Person
	var count = 0
	for _, record := range records {
		markerName := record[0]
		markerValue := record[1]
		if markerValue != "n/a" && markerValue != "" {
			if strings.HasSuffix(markerValue, ".a") ||
				strings.HasSuffix(markerValue, ".g") ||
				strings.HasSuffix(markerValue, ".c") ||
				strings.HasSuffix(markerValue, ".t") {
				markerValue = markerValue[:len(markerValue)-2]
			}
			value, err := strconv.ParseFloat(markerValue, 64)
			if err != nil {
				// This may happen. So just print out a notice.
				fmt.Printf("Error reading YFull marker value in file %s, %s.\n", filename, err)
			} else {
				index, exists := genetic.YFullToIndex(markerName)
				if exists {
					result.YstrMarkers[index] = value
					count++
				} else {
					fmt.Printf("Unknown marker %s found in file %s.\n", markerName, filename)
				}
			}
		}
	}
	// Extract ID and name from filename.
	result.ID = idFromFileName(filepath.Base(filename))
	result.Name = filepath.Base(filename)
	result.Label = stringToLabel(result.ID)
	fmt.Printf("Number of markers for %s: %d\n", result.ID, count)
	return &result, nil
}

// idFromFileName converts the filename of an YFull Y-STR results
// file into the ID of a person.
// A typical filename looks like this: STR_for_YF01234_20160216.csv.
// If the filename is not in the standard format idFromFileName
// tries to extract something useful by stripping the ".csv".
// If this fails the filename itself is returned.
func idFromFileName(filename string) string {
	var result string
	switch {
	case strings.HasPrefix(filename, "STR_for_"):
		end := filename[8:]
		parts := strings.Split(end, "_")
		result = parts[0]
	case strings.HasSuffix(filename, ".csv"):
		result = filename[0 : len(filename)-4]
	default:
		result = filename
	}
	return result
}

// namesWithExt returns the names of all files in a directory
// ending with the extension ext.
// If there are no matching files in the directory
// an empty slice is returned.
func namesWithExt(dirName string, ext string) (filenames []string, err error) {
	filenames = make([]string, 0, 100)
	dir, err := os.Open(dirName)
	if err != nil {
		return filenames, errors.New(fmt.Sprintf("could not open directory, %s\n", err))
	}
	defer dir.Close()

	files, err := dir.Readdirnames(0)
	if err != nil {
		return files, errors.New(fmt.Sprintf("could not read files from directory, %s\n", err))
	}
	for _, filename := range files {
		if filepath.Ext(filename) == ext {
			filenames = append(filenames, filename)
		}
	}
	return filenames, err
}

// WriteDistanceMatrix writes a distance matrix in PHYLIP compatible format.
func WriteDistanceMatrix(filename string, persons []*genetic.Person, matrix *genetic.DistanceMatrix) error {
	// Open file.
	outfile, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer outfile.Close()

	writer := bufio.NewWriter(outfile)
	// Write number of entries
	writer.WriteString(fmt.Sprintf("%d\n", matrix.Size))

	// Write lines
	for row := 0; row < matrix.Size; row++ {
		// Write name
		name := persons[row].Label
		writer.WriteString(name)
		// Write values
		for col := 0; col < matrix.Size; col++ {
			value := strconv.FormatFloat(matrix.Values[row][col], 'f', -1, 64)
			writer.WriteString("\t" + value)
		}
		writer.WriteString("\n")
	}
	err = writer.Flush()
	return err
}

// WritePersonsAsTXT writes person's genetic data to a file.
// The first entry of each line is the person's Label field.
// All entries are separated by tabs so that the content of
// the file can be easily pasted into a spreadsheet.
//
// nMarkers is the number of Y-STR values that is written. This
// is usefull if not all persons have tested for the same number
// of markers.
func WritePersonsAsTXT(filename string, persons []*genetic.Person, nMarkers int) error {
	// Open file.
	outfile, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer outfile.Close()

	writer := bufio.NewWriter(outfile)
	for _, person := range persons {
		writer.WriteString(person.Label)
		for i := 0; i < nMarkers; i++ {
			value := strconv.FormatFloat(person.YstrMarkers[i], 'f', -1, 64)
			writer.WriteString("\t" + value)
		}
		writer.WriteString("\n")
	}
	err = writer.Flush()
	return err
}

// WritePersonsAsHTML writes person's genetic data to a file in HTML format.
//
// nMarkers is the number of Y-STR values that is written. This
// is usefull if not all persons have tested for the same number
// of markers.
func WritePersonsAsHTML(filename string, persons []*genetic.Person, nMarkers int) error {
	modal := persons[0]

	// Open file.
	outfile, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer outfile.Close()

	writer := bufio.NewWriter(outfile)
	// Write header.
	header := "<!DOCTYPE html>\n" +
		"<html lang=\"en\">\n" +
		"<head><title>Y-STR Values</title></head>\n<body>\n<table>\n"
	writer.WriteString(header)

	// Write table data.
	writer.WriteString("<tr>")
	writer.WriteString("<td></td>")
	for i := 0; i < nMarkers; i++ {
		writer.WriteString("<td>" + genetic.YstrMarkerTable[i].InternalName + "</td>")
	}
	writer.WriteString("</tr>\n")
	for _, person := range persons {
		writer.WriteString("<tr>")
		writer.WriteString("<td>" + person.Label + "</td>")
		for i := 0; i < nMarkers; i++ {
			value := strconv.FormatFloat(person.YstrMarkers[i], 'f', -1, 64)
			writer.WriteString("<td style=\"background-color:" +
				colorCode(person.YstrMarkers[i], modal.YstrMarkers[i]) +
				";\">" + value + "</td>")
		}
		writer.WriteString("</tr>\n")
	}
	// Write page end.
	pageEnd := "</table>\n</body>\n</html>"
	writer.WriteString(pageEnd)
	err = writer.Flush()
	return err
}

// colorCode calculates a color for an Y-STR value depending on it's
// distance to the modal value.
// It returns a color string that can be used in CSS stylesheets.
func colorCode(value, modal float64) string {
	if value == 0 {
		return "rgb(242,242,242)"
	}
	colors := [11]string{
		"rgb(0,0,200)",
		"rgb(50,255,255)",
		"rgb(50,255,200)",
		"rgb(50,255,50)",
		"rgb(180,255,180)",
		"rgb(255,255,255)",
		"rgb(255,255,200)",
		"rgb(255,255,100)",
		"rgb(255,200,0)",
		"rgb(255,100,0)",
		"rgb(255,0,0)",
	}
	color := int(value-modal) + 5
	switch {
	case color < 0:
		return colors[0]
	case color >= len(colors):
		return colors[len(colors)-1]
	default:
		return colors[color]
	}
}

// ReadMutationRates reads mutation rates from a file.
// The mutation rates must be provided in JSON format.
func ReadMutationRates(filename string) (genetic.YstrMarkers, error) {
	var result genetic.YstrMarkers
	// Map marker names to indices.
	var names = make(map[string]int)
	for i, _ := range genetic.YstrMarkerTable {
		names[genetic.YstrMarkerTable[i].InternalName] = i
	}
	// Open file.
	infile, err := os.Open(filename)
	if err != nil {
		return result, err
	}
	defer infile.Close()
	// Read JSON from file.
	var untypedJSON interface{}
	decoder := json.NewDecoder(infile)
	err = decoder.Decode(&untypedJSON)
	if err != nil {
		return result, err
	}
	// Use type assertions to get typed data.
	JSON := untypedJSON.(map[string]interface{})
	for key, value := range JSON {
		var markerValue float64
		switch v := value.(type) {
		case int:
			markerValue = float64(v)
		case float64:
			markerValue = v
		default:
			fmt.Printf("ReadMutationRates, unknown value: %v.\n", v)
		}
		// Fill result with values.
		result[names[key]] = markerValue
	}
	return result, nil
}

// WriteMutationRates writes mutation rates to file in JSON format.
func WriteMutationRates(filename string, mutationRates genetic.YstrMarkers) error {
	// Create Json
	var buffer bytes.Buffer
	buffer.WriteString("{")
	n := len(genetic.YstrMarkerTable)
	for i := 0; i < n-1; i++ {
		text := fmt.Sprintf("%q:%G,\n", genetic.YstrMarkerTable[i].InternalName, mutationRates[i])
		buffer.WriteString(text)
	}
	text := fmt.Sprintf("%q:%G}", genetic.YstrMarkerTable[n-1].InternalName, mutationRates[n-1])
	buffer.WriteString(text)

	// Write to file.
	err := ioutil.WriteFile(filename, []byte(buffer.String()), os.ModePerm)
	if err != nil {
		return err
	}
	return nil
}

// stringToLabel transforms a string to a label.
// A label is exactly 10 characters long
// and contains only 8-bit characters.
// All spaces are transformed to underscores.
func stringToLabel(name string) string {
	replacements := map[rune]string{
		' ': "_",
		'Ä': "Ae",
		'ä': "ae",
		'Ü': "Ue",
		'ü': "ue",
		'Ö': "Oe",
		'ö': "oe",
		'ß': "ss",
		'(': "{",
		')': "}",
		':': "_",
		';': "_",
		',': "_",
		'[': "{",
		']': "}",
	}

	// Transform name to a string that contains only non
	// unicode characters and no spaces.
	var buffer bytes.Buffer
	for _, char := range name {
		if newChar, ok := replacements[char]; ok {
			buffer.WriteString(newChar)
		} else if utf8.RuneLen(char) > 1 {
			// Unicode character without replacement. Skip it.
		} else {
			buffer.WriteRune(char)
		}
	}

	// Make sure the new name is 10 characters long.
	newName := buffer.String()
	if len(newName) < 10 {
		// Fill up missing characters with underscores.
		underscores := make([]byte, 10-len(newName))
		for i, _ := range underscores {
			underscores[i] = '_'
		}
		newName = string(underscores) + newName
	} else {
		newName = newName[0:10]
	}
	return newName
}
