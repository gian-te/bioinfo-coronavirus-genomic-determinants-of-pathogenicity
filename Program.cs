using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading;

namespace HighConfidenceAlignmentBlocks
{
    class Program
    {
            
        /*
         * Input:
         * Coded Alignment
         *
         * Output:
         * "High-confidence alignment blocks were identified within the multiple sequence alignment (MSA), which were defined as regions LONGER THAN  15 nt, containing less than 40% gaps IN EACH POSITION"
         */
        static void Main(string[] args)
        {
            List<Tuple<int, int>> positions = new List<Tuple<int, int>>();
            double totalSequencesInHighConfidenceBlocks = 0;
            try
            {
                Dictionary<int, double> gapRatios = new Dictionary<int, double>();

                List<List<string>> genomes = new List<List<string>>();
                List<List<string>> highConfidenceBlocks = new List<List<string>>();

                #region LOAD DATA
                using (var reader = new StreamReader(@"C:\Users\vziex\Desktop\DLSU\BIOINFO\Bioinfo Report 2\new_dataset_final.txt"))
                {

                    while (!reader.EndOfStream)
                    {
                        var line = reader.ReadLine();
                        if (line.StartsWith(">"))
                        {
                            continue;
                        }
                        var values = line.Split(',').ToList();
                        values = values.Skip(1).ToList();
                        genomes.Add(values);
                    }
                }
                #endregion

                // high-confidence alignment blocks
                for (int i = 0; i < 40059;) // for every column
                {
                    var blockStart = i;
                    var column_size = 1; // reset column size
                    i += column_size; // move counter to curent position plus increment size to begin processsing next block after the current block
                    
                 
                    var ratio = ExtractLocalHighConfidenceBlocks(genomes, blockStart, column_size, gapRatios);

                    // the the gap ratio is greater than or equal to 40%, do not add that block
                    if (ratio >= 0.4)
                    {
                        Console.WriteLine("Gap ratio is greater than or equal to 40%, discarding block starting at position {0} to position {1}, Increment Size: {2}", blockStart, blockStart + column_size, column_size);
                   
                    }
                    else
                    {
                        // else if the gap ratio is strictly less than 40%
                        // gap ratio is less than 40% in that column,
                        while (ratio < 0.4)
                        {
                            Console.WriteLine("Continuing...");
                      
                            column_size++; // keep incrementing column size to be added to the block as long as the gap count is less than 40%
                            Console.WriteLine("Processing block starting at position {0} to position {1}, Increment Size: {2}", blockStart, blockStart + column_size, column_size);

                          
                            ratio = ExtractLocalHighConfidenceBlocks(genomes, blockStart, column_size, gapRatios);
                            Console.WriteLine("Gap ratio for column {0} is: {1}% gaps", blockStart + column_size, ratio * 100);


                          
                            Console.WriteLine();
                        }
                        column_size--;
                        if (column_size >= 16) // longer than 15 nucleotides, 1 column of nucleotides = 1 position
                        {
                            // take note of the start and end positions of each high-confidence block
                            positions.Add(new Tuple<int, int>(blockStart, blockStart + column_size));
                            i = blockStart + column_size;
                            blockStart = i;
                            totalSequencesInHighConfidenceBlocks += column_size * 944; // number of sequences = number of columns * 944 rows
                        }

                    }
                }

                Console.WriteLine("Done?");
                Console.WriteLine("Total nucleotides in the blocks: {0}", totalSequencesInHighConfidenceBlocks);

                // should be 53% per paper: "...regions (spanning 53% of the total alignment)..."
                Console.WriteLine("Ratio of sequences in blocks over the entire MSA: {0}%", (totalSequencesInHighConfidenceBlocks / (944 * 40059)) * 100);

                Console.WriteLine("{0} high-confidence blocks detected", positions.Count);
            }
            catch (Exception ex)
            {
                Console.WriteLine("{0} high-confidence blocks detected", positions.Count);

                Console.WriteLine(ex.Message);
            }
        }

        private static double ExtractLocalHighConfidenceBlocks(List<List<string>> genomes, int blockStart, int increment, Dictionary<int, double> gapRatios)
        {
            var genomicRegion = new List<string>();
            double retVal = 0;
            foreach (var genome in genomes) // 944 times
            {
                // assemble the region of a genome from genome[blockStart] until genome[blockStart + increment]
                var region = genome.GetRange(blockStart, increment);

                // join the columns into one substring
                var appendedRegion = string.Join("", region);

                // add the appended region for this genome in a temporary list
                genomicRegion.Add(appendedRegion);

               
            }

            // perform column-wise comparison to check if < 40% gaps IN EACH POSITION (meaning: in each column)

            // the increment value is the number of columns/positions there is in the 944 genomes
            // for each column, 
            for (int columnIndex = 0; columnIndex < increment; columnIndex++)
            {
                var columnStr = "";

                // faster computing, for instance if we are processing 16 columns, no need to re-process the first 15 columns just to process the 16th column
                if (gapRatios.ContainsKey(blockStart + columnIndex))
                {
                    retVal = gapRatios[blockStart + columnIndex];
                    continue;
                }

                // assemble the column at a specific column index by traversing all of the nucleotides in that column.
                foreach (var substring in genomicRegion)
                {
                    // assemble this specific column in this region as a string
                    columnStr += substring[columnIndex].ToString();
                }

                // get the number of sequences in the assembled column
                double totalSequencesInColumn = columnStr.Length;

                // count the number of gaps in the assembled column
                double gapCountInColumn = columnStr.Where(letter => letter.ToString().Equals("0")).Count();

                retVal = gapCountInColumn / totalSequencesInColumn;
                gapRatios.Add(blockStart + columnIndex, retVal);
            }

            return retVal;

        }
    }
}
