using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace CodedAlignmentConverter
{
    class Program
    {
        //identifies which virus type the genomes belongs
        public static string get_virus(string seq)
        {
            if ((seq.Contains("mers")) || (seq.Contains("MERS"))) 
                return "MERS";
            if ((seq.Contains("229e")) || (seq.Contains("229E"))) 
                return "229E";
            if ((seq.Contains("hku1")) || (seq.Contains("HKU1"))) 
                return "HKU1";
            if ((seq.Contains("nl63")) || (seq.Contains("NL63"))) 
                return "NL63";
            if ((seq.Contains("oc43")) || (seq.Contains("OC43"))) 
                return "OC43";
            if ((seq.Contains("sars_cov_2")) || (seq.Contains("SARS-CoV-2"))) 
                return "SARS_2";
            if ((seq.Contains("sars_cov_1")) || (seq.Contains("SARS-CoV-1")) || (seq.Contains("SARS"))) 
                return "SARS_1";
            return "Uncategorized";
        }

        static void Main(string[] args)
        {
            //specify source path here
            string input_path = @"\Dataset\Supp_file1.txt";

            //read the supp_file1.txt contents
            string[] genomes = File.ReadAllLines(input_path);
            List<string> converted = new List<string>();

            //iterates on each line of supp_file1.txt
            foreach (string line in genomes)
            {
                //checks if its the genome name
                if (line.StartsWith(">"))
                {
                    converted.Add(Environment.NewLine);
                    converted.Add(get_virus(line) + ",");
                }
                else
                {
                    string new_value = string.Empty;

                    //iterates on each nucleotide in the line
                    foreach (char chr in line)
                    {
                        //checks if its a gap and replace with 0
                        if (chr == '-')
                            new_value += '0';

                        //check if its a nucleotide and replace with 1
                        else 
                            new_value += '1';
                    }

                    var enumvalue = Enumerable
                        .Range(0, new_value.Length / 1)
                        .Select(x => new_value.Substring(x * 1, 1))
                        .ToList();
                    converted.Add(string.Join(",", enumvalue));
                }
            }

            //specify output path here
            string output_path = @"\bioinfo\supp_file1_bin.txt";

            //write converted binary genome sequence
            File.WriteAllText(output_path, string.Join("", converted));
        }
    }
}
