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

	// Extract persons data from CSV records.
	persons := make([]*genetic.Person, 0, 1000)
	for _, record := range records {
		// Try to extract person data from a record.
		// CSV lines often contain non person data.
		// So an error here happens often and is ignored.
		person, err := personFromFields(record, labelCol)
		if err == nil {
			persons = append(persons, person)
		}
	}
	return persons, nil
}

// personFromFields creates a person from a slice of strings.
//
// The first field must contain the ID of the person. labelIndex
// is the field's index that is used for the person's Label field.
func personFromFields(fields []string, labelIndex int) (*genetic.Person, error) {
	var person genetic.Person
	var err error

	// A valid entry must have at least one Id and 12 Y-STR values.
	// Throw everything else away.
	if len(fields) < 13 {
		return nil, errors.New("not enough data fields")
	}

	// Try to determine the start of the Y-STR values.
	// Results starts with DYS393, which is in the range 9-17.
	resultsStart := 0
	for i, field := range fields {
		dys393, err := strconv.ParseFloat(field, 64)
		if err == nil && dys393 >= 9 && dys393 < 17 {
			resultsStart = i
			break
		}
	}
	// Return if the start of the DYS values could not be found.
	if resultsStart < 1 {
		return nil, errors.New("DYS393 not found")
	}

	// The ID is usually the Kit number.
	person.ID = strings.TrimSpace(fields[0])

	// Return if the ID is too short.
	if len(person.ID) < 1 {
		return nil, errors.New("could not determine person ID")
	}

	person.Label = stringToLabel(strings.TrimSpace(fields[labelIndex]))
	person.YstrMarkers, err = extractYstrMarkers(fields[resultsStart:])
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

	// Make sure that DYS464 contains exactly four values.
	// Many spreadsheets contain four separate columns for DYS464
	// but in Family Tree DNA notation it is denoted like "15-15-17-18".
	DYS464pos := 19
	if strings.Contains(fields[DYS464pos], "-") {
		values := strings.Split(fields[DYS464pos], "-")
		if len(values) > 4 {
			// Copy first four values back into DYS464 field.
			fields[DYS464pos] = strings.Join(values[0:4], "-")
			// Store extra DYS464 values in the upper storage area for DYS464.
			for i := 4; i < len(values); i++ {
				value, err := stringToSTR(values[i])
				if err != nil {
					return markers, err
				}
				markers[genetic.DYS464extStart+i-4] = value
			}
		} else if len(values) < 4 {
			complete := make([]string, 4)
			copy(complete, values)
			fields[DYS464pos] = strings.Join(complete, "-")
		}
	}

	// Concatenate all fields together and separate them by "-", because
	// FamilyTreeDNA uses "-" as a separator for multicopy markers.
	concatFields := strings.Join(fields, "-")

	// Extract the fields again using "-" as a separator.
	// Thus all markers values are separated into individual fields.
	reader := csv.NewReader(strings.NewReader(concatFields))
	reader.Comma = '-'
	strValues, _ := reader.Read()

	// Convert string values to YstrMarkers.
	for i, strValue := range strValues {
		value, err := stringToSTR(strValue)
		if err != nil {
			return markers, err
		}
		markers[i] = value
	}
	return markers, nil
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
	_, err = writer.WriteString(fmt.Sprintf("%d\n", matrix.Size))
	if err != nil {
		return err
	}

	// Write lines
	for row := 0; row < matrix.Size; row++ {
		// Write name
		name := persons[row].Label
		_, err = writer.WriteString(name)
		if err != nil {
			return err
		}
		// Write values
		for col := 0; col < matrix.Size; col++ {
			value := strconv.FormatFloat(matrix.Values[row][col], 'f', -1, 64)
			_, err = writer.WriteString("\t" + value)
			if err != nil {
				return err
			}
		}
		_, err = writer.WriteString("\n")
		if err != nil {
			return err
		}
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
		_, err := writer.WriteString(person.Label)
		for i := 0; i < nMarkers; i++ {
			value := strconv.FormatFloat(person.YstrMarkers[i], 'f', -1, 64)
			_, err = writer.WriteString("\t" + value)
			if err != nil {
				return err
			}
		}
		_, err = writer.WriteString("\n")
		if err != nil {
			return err
		}
	}
	err = writer.Flush()
	return err
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
