// Package phylofriend calculates genetic distances from Y-STR values.
package main

import (
	"flag"
	"fmt"
	"os"
	"strings"

	"code.google.com/p/phylofriend/genetic"
	"code.google.com/p/phylofriend/genfiles"
)

func main() {
	// Command line flags.
	var (
		personsin = flag.String("personsin", "", "Input filename (.txt or .csv).")
		namescol  = flag.Int("namescol", 1, "Column number for names in CSV file.")
		mrin      = flag.String("mrin", "", "Filename for the import of mutation rates.")
		phylipout = flag.String("phylipout", "", "Output filename for PHYLIP distance matrix.")
		txtout    = flag.String("txtout", "", "Output filename for persons in text format.")
		nvalues   = flag.Int("nvalues", 67, "Number of Y-STR values to write into text file.")
		mrout     = flag.String("mrout", "", "Filename for the export of mutation rates.")
		anonymize = flag.Bool("anonymize", false, "Determines if persons should be anonymized.")
		cal       = flag.Float64("cal", 1, "Calibration factor for PHYLIP output.")
		gendist   = flag.Float64("gendist", 25, "Generation distance in years.")
	)
	flag.Parse()

	var (
		persons       []*genetic.Person
		mutationRates genetic.YstrMarkers
		err           error
	)

	// Read mutation rates from file.
	if *mrin != "" {
		mutationRates, err = genfiles.ReadMutationRates(*mrin)
		if err != nil {
			fmt.Printf("Error reading mutation rates %v.\n", err)
			os.Exit(1)
		}
	} else {
		// Use default values.
		mutationRates = genetic.MutationRates
	}

	// Write mutation rates to file.
	if *mrout != "" {
		err = genfiles.WriteMutationRates(*mrout, mutationRates)
		if err != nil {
			fmt.Printf("Error writing mutation rates %v.\n", err)
			os.Exit(1)
		}
	}

	// Read persons from file or exit if no file is provided.
	if *personsin != "" {
		if strings.HasSuffix(strings.ToLower(*personsin), ".csv") {
			persons, err = genfiles.ReadPersonsFromCSV(*personsin, *namescol-1)
		} else {
			persons, err = genfiles.ReadPersonsFromTXT(*personsin)
		}
		if err != nil {
			fmt.Printf("Error loading persons data %v.\n", err)
			os.Exit(1)
		}
	} else {
		// Exit program because all following operations
		// depend on persons data.
		os.Exit(0)
	}

	// Anonymize persons data.
	if *anonymize == true {
		persons = genetic.Anonymize(persons)
	}

	// Write persons data in text format.
	if *txtout != "" {
		err = genfiles.WritePersonsAsTXT(*txtout, persons, *nvalues)
		if err != nil {
			fmt.Printf("Error writing persons data to text file %v.\n", err)
			os.Exit(1)
		}
	}

	// Write distance matrix in phylip compatible format.
	if *phylipout != "" {
		dm := genetic.NewDistanceMatrix(persons, mutationRates, genetic.Distance)
		err = genfiles.WriteDistanceMatrix(*phylipout, persons, dm.Years(*gendist, *cal))
		if err != nil {
			fmt.Printf("Error writing phylip file %v.\n", err)
			os.Exit(1)
		}
	}
}
