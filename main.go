// Package phylofriend calculates genetic distances from Y-STR values.
package main

import (
	"flag"
	"fmt"
	"os"
	"strings"

	"github.com/yogischogi/phylofriend/genetic"
	"github.com/yogischogi/phylofriend/genfiles"
)

func main() {
	// Command line flags.
	var (
		personsin = flag.String("personsin", "", "Input filename (.txt or .csv) or directory.")
		labelcol  = flag.Int("labelcol", 1, "Column number for labels in CSV file.")
		mrin      = flag.String("mrin", "", "Filename for the import of mutation rates.")
		phylipout = flag.String("phylipout", "", "Output filename for PHYLIP distance matrix.")
		txtout    = flag.String("txtout", "", "Output filename for persons in text format.")
		nmarkers  = flag.Int("nmarkers", 0, "Uses only the given number of markers for calculations.")
		mrout     = flag.String("mrout", "", "Filename for the export of mutation rates.")
		anonymize = flag.Bool("anonymize", false, "Anonymizes persons' private data.")
		cal       = flag.Float64("cal", 1, "Calibration factor for PHYLIP output.")
		gentime   = flag.Float64("gentime", 1, "Generation time in years.")
		modal     = flag.Bool("modal", false, "Creates modal haplotype.")
		reduce    = flag.Int("reduce", 1, "Reduces the number of persons (for big trees).")
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
		mutationRates = genetic.DefaultMutationRates()
	}

	// Write mutation rates to file.
	if *mrout != "" {
		err = genfiles.WriteMutationRates(*mrout, mutationRates)
		if err != nil {
			fmt.Printf("Error writing mutation rates %v.\n", err)
			os.Exit(1)
		}
	}

	// Read persons from files or exit if no file or wrong filename is provided.
	if *personsin != "" {
		filenames := strings.Split(*personsin, ",")
		for _, filename := range filenames {
			var pers []*genetic.Person
			fileInfo, err := os.Stat(filename)
			switch {
			case err != nil:
				fmt.Printf("Error, something is wrong with personsin, %v.\n", err)
				os.Exit(1)
			case fileInfo.IsDir():
				pers, err = genfiles.ReadPersonsFromDir(filename)
			case strings.HasSuffix(strings.ToLower(filename), ".csv"):
				pers, err = genfiles.ReadPersonsFromCSV(filename, *labelcol-1)
			default:
				pers, err = genfiles.ReadPersonsFromTXT(filename)
			}
			if err != nil {
				fmt.Printf("Error loading persons data %v.\n", err)
				os.Exit(1)
			}
			persons = append(persons, pers...)
		}
	} else {
		// Exit program because all following operations
		// depend on persons data.
		os.Exit(0)
	}

	// Include only persons who have tested at least for the given number of markers.
	if *nmarkers > 0 {
		persons, err = genetic.ReduceToMarkerSet(persons, *nmarkers)
		if err != nil {
			fmt.Printf("Error reducing persons for the specified number of markers, %v.\n", err)
			os.Exit(1)
		}
	}

	// Reduce amount of data.
	// This is for cases in which the tree gets too large.
	if *reduce > 1 {
		persons, err = genetic.Reduce(persons, *reduce)
		if err != nil {
			fmt.Printf("Error reducing amount of persons, %v.\n", err)
			os.Exit(1)
		}
	}

	// Anonymize persons data.
	if *anonymize == true {
		persons = genetic.Anonymize(persons)
	}

	// Create modal haplotype.
	if *modal == true {
		modal := genetic.ModalHaplotype(persons)
		persons = append(persons, modal)
	}

	// Write persons data in text format.
	if *txtout != "" {
		if *nmarkers > 0 {
			err = genfiles.WritePersonsAsTXT(*txtout, persons, *nmarkers)
		} else {
			err = genfiles.WritePersonsAsTXT(*txtout, persons, genetic.MaxMarkers)
		}
		if err != nil {
			fmt.Printf("Error writing persons data to text file %v.\n", err)
			os.Exit(1)
		}
	}

	// Calculate a distance matrix if the modal value should be
	// calculated or if the matrix should be written to a file.
	var dm *genetic.DistanceMatrix
	if *phylipout != "" || *modal == true {
		dm = genetic.NewDistanceMatrix(persons, mutationRates, genetic.Distance)
		dm = dm.Years(*gentime, *cal)
	}

	// Write distance matrix in phylip compatible format.
	if *phylipout != "" {
		err = genfiles.WriteDistanceMatrix(*phylipout, persons, dm)
		if err != nil {
			fmt.Printf("Error writing PHYLIP file %v.\n", err)
		}
	}

	// Print average distance and standard deviation from modal haplotype.
	if *modal == true {
		// Calculate the average distance from the modal haplotype.
		// The modal haplotype is the last in the distance matrix.
		// The last entry is the distance to itself. So it has to be removed.
		m, s, err := genetic.Average(dm.Values[dm.Size-1][0 : dm.Size-1])
		if err != nil {
			fmt.Printf("Error calculating average and standard deviation, %v.\n", err)
		}
		fmt.Printf("Average distance from modal haplotype: %.2f \u00B1 %.2f\n", m, s)
		fmt.Printf("No correction for Poisson distribution and back mutations.\n")
	}
}
