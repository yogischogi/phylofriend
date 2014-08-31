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
	"strconv"
	"strings"
	"unicode/utf8"

	"code.google.com/p/phylofriend/genetic"
)

// ReadPersonsFromCSV reads persons' data from a CSV file.
// The file may contain comments or other non data rows.
// The function will try to recognize and remove them.
// The first entry of each line containing person data
// must be a unique ID.
// Missing Y-STR values are set to 0.
//
// namesCol is the number of the colum used for the person's
// Name10 field. This will be used as a name when a genetic
// distance matrix is written.
func ReadPersonsFromCSV(filename string, namesCol int) ([]*genetic.Person, error) {
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
		person, err := personFromFields(record, namesCol)
		if err == nil {
			persons = append(persons, person)
		}
	}
	return persons, nil
}

// personFromFields creates a person from a slice of strings.
//
// The first field must contain the ID of the person. nameIndex
// is the field's index that is used for the person's Name10 field.
func personFromFields(fields []string, nameIndex int) (*genetic.Person, error) {
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

	person.Name10 = nameToName10(strings.TrimSpace(fields[nameIndex]))
	person.YstrMarkers, err = extractYstrMarkers(fields[resultsStart:])
	return &person, err
}

// extractYstrMarkers creates YstrMarkers from a slice of text fields.
// The entries must be in FamilyTreeDNA order.
func extractYstrMarkers(fields []string) (genetic.YstrMarkers, error) {
	var markers genetic.YstrMarkers

	// TODO: Make sure that multi copy and palindromic markers
	// have the right length. This might be problematic because
	// some users maintain special spread sheets which hold the
	// data in non FamilyTreeDNA format, e.g. each marker in a
	// single column.

	// Copy all fields together and separate them by -, because
	// FamilyTreeDNA uses - as a separator for multicopy markers.
	var buffer bytes.Buffer
	for i := 0; i < len(fields)-1; i++ {
		buffer.Write([]byte(fields[i]))
		buffer.WriteRune('-')
	}
	buffer.Write([]byte(fields[len(fields)-1]))

	// Extract the fields again using - as a separator.
	// Thus multicopy markers are separated into several fields.
	reader := csv.NewReader(&buffer)
	reader.Comma = '-'
	strValues, _ := reader.Read()

	// Convert string values to YstrMarkers
	for i, strValue := range strValues {
		strValue = strings.TrimSpace(strValue)
		if strValue == "" {
			markers.SetValue(i, 0)
		} else if strValue == "O" {
			markers.SetValue(i, 0)
		} else {
			value, err := strconv.ParseFloat(strValue, 64)
			if err != nil {
				return markers, err
			}
			markers.SetValue(i, value)
		}
	}
	return markers, nil
}

// ReadPersonsFromTXT reads persons' data from a text file.
// The file may contain comment lines starting with //.
// The first entry of each line is used for the person's
// Name10 field that is used for output in distance matrices.
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
		persons[i].Name10 = fields[0]
		nValues := len(fields) - 1
		if genetic.Nmarkers < nValues {
			nValues = genetic.Nmarkers
		}
		for j := 0; j < nValues; j++ {
			value, err := strconv.ParseFloat(fields[j+1], 64)
			if err != nil {
				return persons, err
			}
			persons[i].YstrMarkers.SetValue(j, value)
		}
	}
	return persons, err
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
		name := persons[row].Name10
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
// The first entry of each line is the person's Name10 field.
// All entries are separated by tabs so that the content of
// the file can be easily pasted into a spreadsheet.
//
// nValues is the number of Y-STR values that is written. This
// is usefull if not all persons have tested for the same number
// of markers.
func WritePersonsAsTXT(filename string, persons []*genetic.Person, nValues int) error {
	// Open file.
	outfile, err := os.Create(filename)
	if err != nil {
		return err
	}
	defer outfile.Close()

	writer := bufio.NewWriter(outfile)
	for _, person := range persons {
		_, err := writer.WriteString(person.Name10)
		for i := 0; i < nValues; i++ {
			value := strconv.FormatFloat(person.YstrMarkers.Value(i), 'f', -1, 64)
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
	var mutationRates genetic.YstrMarkers
	infile, err := os.Open(filename)
	if err != nil {
		return mutationRates, err
	}
	defer infile.Close()

	decoder := json.NewDecoder(infile)
	err = decoder.Decode(&mutationRates)
	if err != nil {
		return mutationRates, err
	}
	return mutationRates, nil
}

// WriteMutationRates writes mutation rates to file in JSON format.
func WriteMutationRates(filename string, mutationRates genetic.YstrMarkers) error {
	// Create Json
	buffer := new(bytes.Buffer)
	encoder := json.NewEncoder(buffer)
	encoder.Encode(mutationRates)

	// Convert every entry to a single line for better readabilty.
	outString := strings.Replace(buffer.String(), ",", ",\n", -1)

	// Write to file.
	err := ioutil.WriteFile(filename, []byte(outString), os.ModePerm)
	if err != nil {
		return err
	}
	return nil
}

// NameToName10 transforms a name to a name that is only 10 characters long
// and contains only non unicode characters.
// All spaces are transformed to underscores.
func nameToName10(name string) string {
	replacements := map[rune]string{
		' ': "_",
		'Ä': "Ae",
		'ä': "ae",
		'Ü': "Ue",
		'ü': "ue",
		'Ö': "Oe",
		'ö': "oe",
		'ß': "ss",
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