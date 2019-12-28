using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;

namespace Coil_DHC_Dataset_Maker
{
    public class aaindex
    {
        public static readonly List<aaindex.aaindex_entry> aaindex_entries = aaindex.Load();

        public class aaindex_entry
        {
            public string H_Accession_Number;
            public string D_Data_Description;
            public string R_PMID;
            public string A_Authors;
            public string T_Title_Of_Article;
            public string J_Journal_Reference;
            public List<(string entry_name,double similarity_score)> C_Accession_numbers_of_similar_entries=new List<(string, double)>();
            public List<(char amino_acid,double index_value)> I_Amino_Acid_Index_Data=new List<(char, double)>();
            public List<(char amino_acid,double index_value)> I_Amino_Acid_Index_Data_Normalised=new List<(char, double)>();
        }

        //public static List<(int alphabet_id, string alphabet_name, string alphabet_group, string aaindex_accession_number, descriptive_stats values)> sequence_aaindex_all(string sequence)//, List<aaindex_entry> aaindex_entries = null)
        //{
        //    //

        //    //wrong: if (string.IsNullOrWhiteSpace(sequence)) return new List<descriptive_stats>();

        //    //if (aaindex_entries == null) aaindex_entries = aaindex.Load();

        //    var split_sequence = feature_calcs.split_sequence(sequence);

        //    var sequences = new List<(string name, string sequence)>();
        //    sequences.Add(("unsplit",sequence));
        //    sequences.AddRange(split_sequence.Select(a=>("split",a)).ToList());

        //    var result = new List<(int alphabet_id, string alphabet_name, string alphabet_group, string aaindex_accession_number, descriptive_stats values)>();

        //    var alphabets = feature_calcs.aa_alphabets.ToList();
        //    alphabets.Add((-1,"Overall",new List<string>() { "ARNDCQEGHILKMFPSTWYV" }));
        //    alphabets = alphabets.Where(a => !String.Equals(a.name, "Normal", StringComparison.InvariantCultureIgnoreCase)).ToList();

        //    foreach (var alphabet in alphabets)
        //    {
        //        var result1 = new List<(int alphabet_id, string alphabet_name, string alphabet_group, string aaindex_accession_number, descriptive_stats values)>();

        //        foreach (var group in alphabet.groups)
        //        {
        //            // r = get values for all 500+ aaindex entries

        //            var r = aaindex_entries.SelectMany(a =>
        //            {
        //                //var alphabet_sequence = string.Join("", sequence.Where(b => group.Contains(b)).ToList());
        //                var alphabet_sequences = sequences.Select(b => string.Join("", b.sequence.Where(c => group.Contains(c)).ToList())).ToList();

        //                var result2 = new List<(int alphabet_id, string alphabet_name, string alphabet_group, string aaindex_accession_number, descriptive_stats values)>();
        //                foreach (var alphabet_sequence in alphabet_sequences)
        //                {
        //                    var name = $"{a.H_Accession_Number}_{alphabet.name}_{group}";

        //                    if (string.IsNullOrWhiteSpace(alphabet_sequence))
        //                    {
        //                        var y = descriptive_stats.get_stat_values(null, name);

        //                        var z = (alphabet_id: alphabet.id, alphabet_name: alphabet.name, alphabet_group: group, aaindex_accession_number: a.H_Accession_Number, values: y);

        //                        result2.Add(z);
        //                    }
        //                    else
        //                    {
        //                        var x = sequence_aaindex_entry(aaindex_entries, a.H_Accession_Number, alphabet_sequence);

        //                        var y = descriptive_stats.get_stat_values(x, name);

        //                        var z = (alphabet_id: alphabet.id, alphabet_name: alphabet.name, alphabet_group: group, aaindex_accession_number: a.H_Accession_Number, values: y);

        //                        result2.Add(z);

        //                    }
        //                }

        //                return result2;
        //            }).ToList();

                    
        //            result1.AddRange(r);
        //        }

        //        result.AddRange(result1);
        //    }

        //    return result;

        //}

        public static List<(char amino_acid, double value)> sequence_aaindex_entry(/*List<aaindex_entry> aaindex_list, */string accession_number, string sequence)
        {
            // returns a list of values for each sequence

            if (string.IsNullOrWhiteSpace(sequence)) return new List<(char amino_acid, double value)>();

            //if (aaindex_list == null && aaindex.aaindex_entries != null) aaindex_list = aaindex.aaindex_entries;

            var aaindex_entry = aaindex_entries.FirstOrDefault(a => a.H_Accession_Number == accession_number);
            
            var result = sequence.Select(aa => (amino_acid:aa, value:aaindex_entry.I_Amino_Acid_Index_Data.First(b => b.amino_acid == aa).index_value)).ToList();

            return result;
        }
        

        public static List<aaindex_entry> Load(string aaindex1_file= @"C:\betastrands_dataset\aaindex\aaindex1.txt", bool remove_incomplete_entries = false)
        {
            var result = new List<aaindex_entry>();
            var lines = File.ReadAllLines(aaindex1_file);
            aaindex_entry entry=null;
            
            var i_amino_acid_order = "";
            var lastcode = "";
            for (var lines_index = 0; lines_index < lines.Length; lines_index++)
            {
                var line = lines[lines_index];
                if (line.Length == 0) continue;
                var code = char.IsWhiteSpace(line, 0) ? null : line.Split().FirstOrDefault();
                if (string.IsNullOrWhiteSpace(code)) code = lastcode;
                lastcode = code;
                if (line.Length == 1) continue;
                if (code == "//" || entry == null || lines_index == 0)
                {
                    if (lines_index >= lines.Length - 2) break;
                    entry = new aaindex_entry();
                    result.Add(entry);
                    i_amino_acid_order = "";
                    if (code == "//")
                    {
                        lastcode = "";
                        code = "";
                    }
                }

                var line_data = line.Substring(2).Trim();


                switch (code)
                {
                    case "H":
                        if (!string.IsNullOrEmpty(entry.H_Accession_Number)) entry.H_Accession_Number += Environment.NewLine;
                        entry.H_Accession_Number += line_data;
                        break;
                    case "D":

                        if (!string.IsNullOrEmpty(entry.D_Data_Description)) entry.D_Data_Description += Environment.NewLine;
                        entry.D_Data_Description += line_data;
                        break;
                    case "R":

                        if (!string.IsNullOrEmpty(entry.R_PMID)) entry.R_PMID += Environment.NewLine;
                        entry.R_PMID += line_data;
                        break;
                    case "A":

                        if (!string.IsNullOrEmpty(entry.A_Authors)) entry.A_Authors += Environment.NewLine;
                        entry.A_Authors += line_data;
                        break;
                    case "T":

                        if (!string.IsNullOrEmpty(entry.T_Title_Of_Article)) entry.T_Title_Of_Article += Environment.NewLine;
                        entry.T_Title_Of_Article += line_data;
                        break;
                    case "J":

                        if (!string.IsNullOrEmpty(entry.J_Journal_Reference)) entry.J_Journal_Reference += Environment.NewLine;
                        entry.J_Journal_Reference += line_data;
                        break;

                    case "*":
                        break;
                    case "C":
                        var c = line_data.Split(new char[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries).ToList();
                        while (c.Count > 0)
                        {
                            entry.C_Accession_numbers_of_similar_entries.Add((c[0], double.Parse(c[1])));
                            c = c.Skip(2).ToList();
                        }

                        break;
                    case "I":
                        if (string.IsNullOrWhiteSpace(i_amino_acid_order))
                        {
                            //if (entry.H_Accession_Number == "OOBM770105")
                            //{
                            //    Program.WriteLine();
                            //}
                            i_amino_acid_order = string.Join("", line_data.Split(new char[] {' ', '\t', '/', '\\'}, StringSplitOptions.RemoveEmptyEntries).Where((a, i) => i % 2 == 0).ToList()) + string.Join("", line_data.Split(new char[] {' ', '\t', '/', '\\'}, StringSplitOptions.RemoveEmptyEntries).Where((a, i) => i % 2 != 0).ToList());
                            //A/L     R/K     N/M     D/F     C/P     Q/S     E/T     G/W     H/Y     I/V   
                        }
                        else
                        {
                            var nums = line_data.Split(new char[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries).Select(a => a.Trim()).ToList();
                            foreach (var num in nums)
                            {
                                if (double.TryParse(num, out var numd))
                                {
                                    entry.I_Amino_Acid_Index_Data.Add((i_amino_acid_order[0], numd));
                                }
                                else if (!remove_incomplete_entries)
                                {
                                    entry.I_Amino_Acid_Index_Data.Add((i_amino_acid_order[0], 0));
                                }

                                i_amino_acid_order = i_amino_acid_order.Substring(1);
                            }
                        }

                        break;
                }
            }

            if (remove_incomplete_entries) result = result.Where(a => a.I_Amino_Acid_Index_Data.Count >= 20).ToList();

            var normalise = true;

            if (normalise)
            {

                // normalise
                foreach (var r in result)
                {
                    var min = r.I_Amino_Acid_Index_Data.Min(a => a.index_value);
                    var max = r.I_Amino_Acid_Index_Data.Max(a => a.index_value);

                    r.I_Amino_Acid_Index_Data_Normalised = r.I_Amino_Acid_Index_Data.Select(a => (a.amino_acid, scale_value(a.index_value, min, max))).ToList();
                }

                // make overall average from normalised data
                var x = result.SelectMany(a => a.I_Amino_Acid_Index_Data_Normalised).GroupBy(a => a.amino_acid).Select(a => (amino_acid:a.Key, index_value:a.Select(b=>b.index_value).Average())).ToList();

                var aaindex_overall_average_str = "aaindex_overall_average";
                var aaindex_overall_average = new aaindex_entry()
                {
                    A_Authors = aaindex_overall_average_str,
                    H_Accession_Number = aaindex_overall_average_str,
                    C_Accession_numbers_of_similar_entries = new List<(string entry_name, double similarity_score)>(),
                    D_Data_Description = aaindex_overall_average_str,
                    I_Amino_Acid_Index_Data = x,
                    I_Amino_Acid_Index_Data_Normalised = x,
                    J_Journal_Reference = aaindex_overall_average_str,
                    R_PMID = aaindex_overall_average_str,
                    T_Title_Of_Article = aaindex_overall_average_str
                };
       
                result.Add(aaindex_overall_average);
            }

            return result;
        }

        public static double scale_value(double value, double min, double max)
        {
            var x = value - min;
            var y = max - min;

            if (x == 0) return 0;
            if (y == 0) return value;

            var scaled = x / y;

            return scaled;
        }
    }
}
