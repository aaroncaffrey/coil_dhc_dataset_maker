using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace Coil_DHC_Dataset_Maker
{
    public partial class feature_calcs
    {
        //public static List<Tuple<char, char, int>> SequencePropensityPairs(string sequence)
        //{
        //    var sequencePad = sequence + "__";
        //    var subsequences = sequence.Select((a, i) => sequencePad.Substring(i, 2)).ToList();
        //    subsequences.AddRange(sequence.Select((a, i) => sequencePad.Substring(i + 1, 2)).ToList());

        //    var aa = new char[] { 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' };

        //    var scores = new List<Tuple<char, char, int>>();
        //    foreach (var a1 in aa)
        //    {
        //        foreach (var a2 in aa)
        //        {
        //            var c = subsequences.Count(b => b == "" + a1 + a2);

        //            scores.Add(new Tuple<char, char, int>(a1, a2, c));
        //        }
        //    }

        //    return scores;
        //}

        //public static List<Tuple<char, int>> SequencePropensity(string sequence)
        //{
        //    var aa = new char[] { 'A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V' };

        //    var prop = aa.Select(a => new Tuple<char, int>(a, sequence.Count(b => b == a))).ToList();

        //    return prop;
        //}

        public static double transform_value(double value, double length, bool sqrt, bool as_dist)
        {
            if (value == 0.0)
            {
                return value;
            }

            if (as_dist && length != 0)
            {
                value = value / length;
            }

            if (sqrt)
            {
                value = (double)Math.Sqrt((double)value);
            }

            return value;
        }

        //public static List<(char ss_type, double value)> ss_distribution(string ss_sequence, bool sqrt, bool as_dist, ss_types ss_type)// = ss_types.DSSP3)
        //{
        //    const string dssp3_types = "EHC";
        //    const string dssp_types = "HBEGITSC";

        //    const string stride3_types = "EHC";
        //    const string stride_types = "HBEGITbC";

        //    var ss_codes = "";

        //    if (ss_type == ss_types.DSSP)
        //    {
        //        ss_codes = dssp_types;
        //    }
        //    else if (ss_type == ss_types.DSSP3)
        //    {
        //        ss_codes = dssp3_types;
        //        ss_sequence = Atom.secondary_structure_state_reduction(ss_sequence);
        //    }
        //    else if (ss_type == ss_types.STRIDE)
        //    {
        //        ss_codes = stride_types;
        //    }
        //    else if (ss_type == ss_types.STRIDE3)
        //    {
        //        ss_codes = stride3_types;
        //        ss_sequence = Atom.secondary_structure_state_reduction(ss_sequence);
        //    }


        //    var d = ss_codes.Select(a => (a, string.IsNullOrWhiteSpace(ss_sequence) ? (double)0 : transform_value(ss_sequence.Count(b => a == b), ss_sequence.Length, sqrt, as_dist))).ToList();

        //    return d;
        //}



        //public static List<(char aa_type, double value)> aa_distribution(string sequence, bool sqrt, bool as_dist)
        //{
        //    // 20 features
        //    var prop = aa_list.Select(aa => (aa_type: aa, value: string.IsNullOrWhiteSpace(sequence) ? (double)0 :
        //        transform_value(sequence.Count(seq_aa => seq_aa == aa), sequence.Length, sqrt, as_dist))).ToList();
        //    return prop;
        //}

        //public static List<(string pp_group_name, double value)> aa_group_distribution(string sequence, bool sqrt, bool as_dist)
        //{
        //    // 13 features
        //    var prop = aa_chem_groups.Select(group =>
        //        (pp_group_name: group.group_name, value: string.IsNullOrWhiteSpace(sequence) ? (double)0 :
        //            transform_value(sequence.Count(seq_aa => group.aa_list.Contains(seq_aa)), sequence.Length, sqrt, as_dist))).ToList();
        //    return prop;
        //}

        //public List<(int id, string name, double[] values)> feature_AAComposition(string seq, bool sqrt, bool as_dist)
        //{

        //    var result = new List<(int id, string name, double[] values)>();
        //    for (var alphabet_index = 0; alphabet_index < alphabets.Count; alphabet_index++)
        //    {
        //        var alphabet = alphabets[alphabet_index];
        //        var alphabet_result = new double[alphabet.groups.Count];

        //        for (var group_index = 0; group_index < alphabet.groups.Count; group_index++)
        //        {
        //            var group = alphabet.groups[group_index];

        //            var value = (double) seq.Count(a => group.Contains(a));
        //            value = transform_value(value, seq.Length, sqrt, as_dist);
        //            alphabet_result[group_index] = value;
        //        }

        //        result.Add((alphabet.id,alphabet.name,alphabet_result));
        //    }

        //    return result;
        //}

        public static readonly List<(int id, string name, List<(string group_name, string group_amino_acids)> groups)> aa_alphabets = new List<(int id, string name, List<(string group_name, string group_amino_acids)> groups)>()
        {
            //(-1, "Overall",new List<string>(){
            //    "ARNDCQEGHILKMFPSTWYV"
            //}),
            (0, "Normal",new List<(string group_name, string group_amino_acids)>(){
                ("A","A"), ("R","R"), ("N","N"), ("D","D"), ("C","C"), ("Q","Q"), ("E","E"),
                ("G","G"), ("H","H"), ("I","I"), ("L","L"), ("K","K"), ("M","M"), ("F","F"),
                ("P","P"), ("S","S"), ("T","T"), ("W","W"), ("Y","Y"), ("V","V")
                // invalid: B J O U X Z
            }),
            (1, "Physicochemical",new List<(string group_name, string group_amino_acids)>(){
                ("HydrophobicSmall_AVFPMILW","AVFPMILW"),
                ("Negative_DE","DE"),
                ("Positive_RK","RK"),
                ("HydroxylAmineBasic_STYHCNGQ","STYHCNGQ")
            }),
            (2, "Hydrophobicity",new List<(string group_name, string group_amino_acids)>(){
                ("Alaphatic_LAGVIP","LAGVIP"),
                ("Aromatic_FYW","FYW"),
                ("Polar_DENQRHSTK","DENQRHSTK"),
                ("Sulphuric_CM","CM")
            }),
            (3, "UniProtKb",new List<(string group_name, string group_amino_acids)>(){
                ("LAGVIP","LAGVIP"),
                ("Negative_DE","DE"),
                ("Hydroxylic_ST","ST"),
                ("Positive_HKR","HKR"),
                ("Aromatic_FYW","FYW"),
                ("Acidic_NQ","NQ"),
                ("Sulphuric_CM", "CM")
            }),
            (4, "PdbSum",new List<(string group_name, string group_amino_acids)>(){
                ("Positive_HKR","HKR"),
                ("Negative_DE","DE"),
                ("Polar_STNQ","STNQ"),
                ("Hydrophobic_AVLIM","AVLIM"),
                ("Aromatic_FYW","FYW"),
                ("Small_PG","PG"),
                ("Sulphuric_C","C")
            }),
            (5, "Venn1",new List<(string group_name, string group_amino_acids)>(){
                ("Aliphatic_ILV","ILV"),
                ("Hydroxylic_TS","TS"),
                ("Tiny_AGCS","AGCS"),
                ("Small_VPAGCSTDN","VPAGCSTDN"),
                ("Acidic_NQ","NQ"),
                ("Positive_HKR","HKR"),
                ("Negative_DE","DE"),
                ("Polar_CSTNQDEHKRYW","CSTNQDEHKRYW"),
                ("Nonpolar_ILVAMFGP","ILVAMFGP"),
                ("Charged_DEHKR","DEHKR"),
                ("Hydrophobic_ILVACTMFYWHK","ILVACTMFYWHK"),
                ("Aromatic_FYWH","FYWH"),
                ("Sulphuric_MC","MC")
            }),
            (6, "Venn2",new List<(string group_name, string group_amino_acids)>(){
                ("Tiny_ACGST","ACGST"),
                ("Small_ACDGNPSTV","ACDGNPSTV"),
                ("Aliphatic_AILV","AILV"),
                ("Aromatic_FHWY","FHWY"),
                ("Nonpolar_ACFGILMPVWY","ACFGILMPVWY"),
                ("Polar_DEHKNQRST","DEHKNQRST"),
                ("Charged_DEHKR","DEHKR"),
                ("Positive_HKR","HKR"),
                ("Negative_DE","DE")
            }),
            (7, "TableVenn",new List<(string group_name, string group_amino_acids)>(){
                ("Acyclic_ARNDCEQGILKMSTV","ARNDCEQGILKMSTV"),
                ("Aliphatic_AGILV","AGILV"),
                ("Aromatic_HFWY","HFWY"),
                ("Buried_ACILMFWV","ACILMFWV"),
                ("Charged_RDEHK","RDEHK"),
                ("Cyclic_HFPWY","HFPWY"),
                ("Hydrophobic_AGILMFPWYV","AGILMFPWYV"),
                ("Large_REQHILKMFWY","REQHILKMFWY"),
                ("Medium_NDCPTV","NDCPTV"),
                ("Negative_DE","DE"),
                ("Neutral_ANCQGHILMFPSTWYV","ANCQGHILMFPSTWYV"),
                ("Polar_RNDCEQHKST","RNDCEQHKST"),
                ("Positive_RHK","RHK"),
                ("Small_AGS","AGS"),
                ("Surface_RNDEQGHKPSTY","RNDEQGHKPSTY"),
            })
        };

        public static readonly List<(int id, string name, List<(string group_name, string group_amino_acids)> groups)> aa_alphabets_foldx =
            aa_alphabets.Select(a =>
            {
                var x = (a.id, a.name, a.groups.Select(b =>
                {
                    var y = b.group_amino_acids + string.Join("", foldx_caller.foldx_residues_aa_mutable.Where(d => d.standard_aa_code1 != d.foldx_aa_code1 && b.group_amino_acids.Contains(d.standard_aa_code1)).ToList());
                    return (b.group_name,y);
                }).ToList());
                return x;
            }).ToList();

        public static readonly List<(int id, string name, List<(string group_name, string group_amino_acids)> groups)> aa_alphabets_inc_overall = new List<(int id, string name, List<(string group_name, string group_amino_acids)> groups)>(aa_alphabets)
        {
            (-1, "Overall",new List<(string group_name, string group_amino_acids)> (){
                ("Overall_ARNDCQEGHILKMFPSTWYV","ARNDCQEGHILKMFPSTWYV")
            }),
        };


        public static readonly List<(int id, string name, List<(string group_name, string group_amino_acids)> groups)> aa_alphabets_inc_overall_foldx =
            aa_alphabets_inc_overall.Select(a =>
            {
                var x = (a.id, a.name, a.groups.Select(b =>
                {
                    var y = b.group_amino_acids + string.Join("", foldx_caller.foldx_residues_aa_mutable.Where(d => d.standard_aa_code1 != d.foldx_aa_code1 && b.group_amino_acids.Contains(d.standard_aa_code1)).ToList());
                    return (b.group_name,y);
                }).ToList());
                return x;
            }).ToList();

        public static readonly List<(int id, string name, List<(string group_name, string group_amino_acids)> groups)> ss_alphabets = new List<(int id, string name, List<(string group_name, string group_amino_acids)> groups)>()
        {
            (-1, "SS_HEC"/*T*/,new List<(string group_name, string group_amino_acids)>(){
                ("Helix_H","H"),
                ("Strand_E","E"),
                ("Coil_C","C"),
                //("Turn_T","T")
            }),
        };

        //        (0, "Aliphatic", "ILV"),
        //        (1, "Hydroxylic", "TS"),
        //        (2, "Tiny", "AGCS"),
        //        (3, "Small", "VPAGCSTDN"),
        //        (4, "Acidic", "NQ"),
        //        (5, "Positive", "HKR"),
        //        (6, "Negative", "DE"),
        //        (7, "Polar", "CSTNQDEHKRYW"),
        //        (8, "Non-polar", "ILVAMFGP"),
        //        (9, "Charged", "DEHKR"),
        //        (10, "Hydrophobic", "ILVACTMFYWHK"),
        //        (11, "Aromatic", "FYWH"),
        //        (12, "Sulphuric", "MC"),

        public class named_double
        {
            public string name;
            public double value;
        }

        public static string[] split_sequence(string seq, int sections = 3, bool distribute = false)
        {
            if (seq == null) seq = "";

            var x = split_sequence(seq.ToList(), sections, 0, distribute);
            var y = x.Select(a => string.Join("", a)).ToArray();
            return y;
        }


        public static List<List<T>> split_sequence<T>(List<T> seq, int sections, int divisible = 0, bool distribute = false)
        {
            var r = new List<List<T>>();
            for (var i = 0; i < sections; i++)
            {
                r.Add(new List<T>());
            }

            if (seq == null || seq.Count == 0)
            {
                return r;
            }

            var len = seq.Count;

            if (sections > len) sections = len;

            var div = len / sections;
            var mod = len % sections;
            var mod_start_index = ((sections / 2) - (mod / 2));
            var mod_end_index = mod_start_index + (mod - 1);
            var centre = (sections / 2);
            var taken = 0;
            for (var i = 0; i < sections; i++)
            {
                if (taken >= len) break;
                var take = div + (distribute && i >= mod_start_index && i <= mod_end_index ? 1 : 0) + (!distribute && i == centre ? mod : 0);

                if (divisible != 0)
                {
                    var rem = take % divisible;

                    take += rem;
                }

                var x = seq.Skip(taken).Take(take).ToList();

                r[i].AddRange(x);
                taken += x.Count;
            }

            return r;
        }


        public enum seq_type
        {
            amino_acid_sequence,
            secondary_structure_sequence
        }

        public static 
            List<( 
                int alphabet_id, 
                string alphabet_name, 
                named_double[][] motifs, // A B C ... AA AB AC ... AAA AAB AAC ...
                named_double[][] motifs_binary, 
                //named_double[] saac, // saac
                //named_double[] saac_binary, 
                named_double[] oaac, // oaac
                named_double[] oaac_binary, 
                named_double[] average_seq_positions, // average position of each AA in sequence,
                named_double[][] dipeptides, // AA A?A A??A AB A?B A??B // the average distance between A and B.  What about between AX and BX...?  800 features with full alphabet, or 32 features with 4 letter alphabet... but with only 1 distance, or 64 with 2 distance
                named_double[][] dipeptides_binary, // AA A?A A??A AB A?B A??B
                named_double[] average_dipeptide_distance // average distance between each AA
            )> 
            feature_pse_aac(string seq, seq_type seq_type, pse_aac_options pse_aac_options, bool sqrt, bool as_dist)
        {
            // number of distinct amino acids in sequence
            // number of continuous group amino acids

            // Sequence features:
            // 1. Amino Acid Composition
            // 2. Amino Acid Order Context (How many times do A and B appear at X intervals?)
            // 3. Amino Acid Order Distance (What is the average distance between A and B?)
            // 4. Amino Acid Average Position (What is the average % position of each amino acid within the sequence)?
            // 5. motifs of length 1, 2 or 3
            var alphabets = seq_type == seq_type.amino_acid_sequence ? aa_alphabets : ss_alphabets;

            //var seq_split = split_sequence(seq);

            var indexes0based = alphabets.Select(alphabet => alphabet.groups.Select(alphabet_group => (alphabet_group: alphabet_group, indexofall: seq.IndexOfAll(alphabet_group.group_amino_acids.ToCharArray()).ToList())).ToList()).ToList();

            //var split_indexes0based = alphabets.Select(alphabet => alphabet.groups.Select(alphabet_group => (alphabet: alphabet, alphabet_group: alphabet_group, indexofall: seq_split.Select(a_seq_split => (alphabet: alphabet, alphabet_group: alphabet_group, a_seq_split: a_seq_split, indexofall: a_seq_split.IndexOfAll(alphabet_group.group_amino_acids.ToCharArray()).ToList())).ToList())).ToList()).ToList();


            var a_tasks = new List<Task<List<(
                int alphabet_id,
                string alphabet_name, 
                named_double[][] motifs, // A B C ... AA AB AC ... AAA AAB AAC ...
                named_double[][] motifs_binary,
                //named_double[] saac, // saac
                //named_double[] saac_binary,
                named_double[] oaac, // oaac
                named_double[] oaac_binary, 
                named_double[] average_seq_positions, // average position of each AA in sequence,
                named_double[][] dipeptides, // AA A?A A??A AB A?B A??B
                named_double[][] dipeptides_binary, // AA A?A A??A AB A?B A??B
                named_double[] average_dipeptide_distance // average distance between each AA
                )>>>();

            //var alpha_tasks = new List<Task>();
            var sequence_relation_levels = 3;
            for (var ext_alphabet_index = 0; ext_alphabet_index < alphabets.Count; ext_alphabet_index++)
            {
                var alphabet_index = ext_alphabet_index;

                var a_task = Task.Run(() =>
                {
                    var result = new List<(
                        int alphabet_id, 
                        string alphabet_name, 
                        named_double[][] motifs, // A B C ... AA AB AC ... AAA AAB AAC ...
                        named_double[][] motifs_binary,
                        //named_double[] saac, // saac
                        //named_double[] saac_binary,
                        named_double[] oaac, // oaac
                        named_double[] oaac_binary,
                        named_double[] average_seq_positions, // average position of each AA in sequence,
                        named_double[][] dipeptides, // AA A?A A??A AB A?B A??B
                        named_double[][] dipeptides_binary, // AA A?A A??A AB A?B A??B
                        named_double[] average_dipeptide_distance // average distance between each AA
                        )>();

                    var tasks = new List<Task>();
                    var alphabet = alphabets[alphabet_index];

                    var alphabet_indexes0based = indexes0based[alphabet_index];
                    //var split_alphabet_indexes0based = split_indexes0based[alphabet_index];







                    named_double[][] motifs = null;
                    if (pse_aac_options.motifs || pse_aac_options.motifs_binary)
                    {
                        

                        var task = Task.Run(() =>
                        {
                            motifs = new named_double[sequence_relation_levels][];

                            for (var i = 0; i < sequence_relation_levels; i++)
                            {
                                motifs[i] = new named_double[(int) Math.Pow(alphabet.groups.Count, i + 1)];
                            }

                            var motif_index1 = -1;
                            var motif_index2 = -1;
                            var motif_index3 = -1;
                            for (var alphabet_group_index1 = 0; alphabet_group_index1 < alphabet.groups.Count; alphabet_group_index1++)
                            {
                                motif_index1++;
                                var alphabet_group1 = alphabet.groups[alphabet_group_index1];
                                var motif_expressions1 = alphabet_group1.group_amino_acids.Select(a => a).ToList();
                                var count1 = string.IsNullOrEmpty(seq) ? 0d : (double) seq.Count(a => motif_expressions1.Contains(a));
                                count1 = transform_value(count1, seq == null ? 0 : seq.Length, sqrt, as_dist);
                                motifs[0][motif_index1] = new named_double() {name = alphabet_group1.group_amino_acids, value = count1};

                                // calc motif len 1 
                                for (var alphabet_group_index2 = 0; alphabet_group_index2 < alphabet.groups.Count; alphabet_group_index2++)
                                {
                                    motif_index2++;
                                    var alphabet_group2 = alphabet.groups[alphabet_group_index2];
                                    var motif_expressions2 = alphabet_group1.group_amino_acids.SelectMany(a => alphabet_group2.group_amino_acids.Select(b => $"{a}{b}").ToList()).ToList();
                                    var count2 = string.IsNullOrEmpty(seq) ? 0d : (double) motif_expressions2.Sum(m => seq.IndexOfAll(m).Count);
                                    count2 = transform_value(count2, seq == null ? 0 : seq.Length, sqrt, as_dist);
                                    motifs[1][motif_index2] = new named_double() {name = $"{alphabet_group1.group_amino_acids}_{alphabet_group2.group_amino_acids}", value = count2};

                                    // calc motif len 2

                                    for (var alphabet_group_index3 = 0; alphabet_group_index3 < alphabet.groups.Count; alphabet_group_index3++)
                                    {
                                        motif_index3++;
                                        var alphabet_group3 = alphabet.groups[alphabet_group_index3];

                                        var motif_expressions3 = alphabet_group1.group_amino_acids.SelectMany(a => alphabet_group2.group_amino_acids.SelectMany(b => alphabet_group3.group_amino_acids.Select(c => $"{a}{b}{c}").ToList()).ToList()).ToList();

                                        var count3 = string.IsNullOrEmpty(seq) ? 0d : (double) motif_expressions3.Sum(m => seq.IndexOfAll(m).Count);
                                        count3 = transform_value(count3, seq == null ? 0 : seq.Length, sqrt, as_dist);

                                        motifs[2][motif_index3] = new named_double() {name = $"{alphabet_group1.group_amino_acids}_{alphabet_group2.group_amino_acids}_{alphabet_group3.group_amino_acids}", value = count3};

                                        // calc motif len 3
                                    }
                                }
                            }
                        });
                        tasks.Add(task);
                    }


                    named_double[] oaac = null;
                    named_double[] average_seq_positions = null;
                    //named_double[] saac = null;
                    //named_double[][] saac_jagged = null;
                    if (pse_aac_options.oaac || pse_aac_options.oaac_binary || /*pse_aac_options.saac || pse_aac_options.saac_binary ||*/ pse_aac_options.average_seq_position)
                    {
                        
                        

                        var task = Task.Run(() =>
                        {
                            if (pse_aac_options.oaac || pse_aac_options.oaac_binary)
                            {
                                oaac = new named_double[alphabet.groups.Count];
                            }

                            if (pse_aac_options.average_seq_position)
                            {
                                average_seq_positions = new named_double[alphabet.groups.Count];
                            }

                            //if (pse_aac_options.saac || pse_aac_options.saac_binary)
                            //{
                            //    saac_jagged = new named_double[alphabet.groups.Count][];
                            //}

                            for (var alphabet_group_index = 0; alphabet_group_index < alphabet.groups.Count; alphabet_group_index++)
                            {
                                var alphabet_group = alphabet.groups[alphabet_group_index];

                                var alphabet_group_indexes0based = alphabet_indexes0based[alphabet_group_index].indexofall;

                                if (pse_aac_options.oaac || pse_aac_options.oaac_binary)
                                {
                                    var group_instance_count = (double) alphabet_group_indexes0based.Count;
                                    
                                    group_instance_count = transform_value(group_instance_count, string.IsNullOrEmpty(seq) ? 0 : seq.Length, sqrt, as_dist);

                                    oaac[alphabet_group_index] = new named_double() {name = alphabet_group.group_amino_acids, value = group_instance_count};
                                }

                                //if (pse_aac_options.saac || pse_aac_options.saac_binary)
                                //{

                                //    var split_alphabet_group_indexes0based = split_alphabet_indexes0based[alphabet_group_index];
                                    
                                //    var split_group_instance_count = split_alphabet_group_indexes0based.indexofall.Select(a => (double) a.indexofall.Count).ToList();

                                //    split_group_instance_count = split_group_instance_count.Select((a,i) => transform_value(a, string.IsNullOrEmpty(seq) ? 0 : seq_split[i].Length, sqrt, as_dist)).ToList();

                                //    saac_jagged[alphabet_group_index] = split_group_instance_count.Select((a, i) => new named_double() {name = $"{alphabet_group.group_amino_acids}_{i}", value = a}).ToArray();
                                //}

                                if (pse_aac_options.average_seq_position)
                                {
                                    var average_position = string.IsNullOrEmpty(seq) ? 0 : (alphabet_group_indexes0based.Count > 0 ? (alphabet_group_indexes0based.Average() / (seq.Length - 1)) : 0.5);

                                    average_seq_positions[alphabet_group_index] = new named_double() {name = alphabet_group.group_amino_acids, value = average_position};
                                }
                            }


                            //if (pse_aac_options.saac || pse_aac_options.saac_binary)
                            //{
                            //    saac = saac_jagged.SelectMany(a => a).ToArray();
                            //}
                        });
                        tasks.Add(task);
                    }



                    named_double[] average_dipeptide_distance = null;
                    named_double[][] dipeptides = null;
                    if (pse_aac_options.average_dipeptide_distance || pse_aac_options.dipeptides || pse_aac_options.dipeptides_binary)
                    {
                        

                        var task = Task.Run(() =>
                        {
                            var alphabet_aa_pairs = alphabet.groups.SelectMany((a, x) => alphabet.groups.Select((b, y) => (x, y)).ToList()).Where((c, y) => c.x <= c.y).ToList();


                            if (pse_aac_options.average_dipeptide_distance)
                            {
                                average_dipeptide_distance = new named_double[alphabet_aa_pairs.Count];
                            }



                            if (pse_aac_options.dipeptides || pse_aac_options.dipeptides_binary)
                            {
                                dipeptides = new named_double[sequence_relation_levels][];
                                for (var i = 0; i < dipeptides.Length; i++)
                                {
                                    dipeptides[i] = new named_double[alphabet_aa_pairs.Count];
                                }
                            }



                            for (var alphabet_pairs_index = 0; alphabet_pairs_index < alphabet_aa_pairs.Count; alphabet_pairs_index++)
                            {
                                var alphabet_pair = alphabet_aa_pairs[alphabet_pairs_index];

                                var group1 = alphabet.groups[alphabet_pair.x];
                                var group2 = alphabet.groups[alphabet_pair.y];

                                var indexes1 = alphabet_indexes0based[alphabet_pair.x].indexofall;
                                var indexes2 = alphabet_indexes0based[alphabet_pair.y].indexofall;

                                var distances = indexes1.SelectMany(a => indexes2.Select(b => Math.Abs(a - b)).ToList()).ToList();

                                if (pse_aac_options.average_dipeptide_distance)
                                {
                                    var distance_value = distances.Count > 0 ? distances.Average() : 0;
                                    distance_value = transform_value(distance_value, string.IsNullOrEmpty(seq) ? 0 : seq.Length, sqrt, as_dist);
                                    average_dipeptide_distance[alphabet_pairs_index] = new named_double() {name = $"{group1.group_amino_acids}_{group2.group_amino_acids}", value = distance_value};
                                }

                                if (pse_aac_options.dipeptides || pse_aac_options.dipeptides_binary)
                                {
                                    for (var n_index = 0; n_index < sequence_relation_levels; n_index++)
                                    {
                                        var n = n_index + 1;
                                        var context_n = (double) distances.Count(a => a == n);
                                        context_n = transform_value(context_n, string.IsNullOrEmpty(seq) ? 0 : seq.Length, sqrt, as_dist);
                                        dipeptides[n_index][alphabet_pairs_index] = new named_double() {name = $"{group1.group_amino_acids}_{group2.group_amino_acids}", value = context_n};
                                    }
                                }
                            }
                        });
                        tasks.Add(task);
                    }

                    Task.WaitAll(tasks.ToArray<Task>());


                    
                    
                    
                    

                    var tasks2 = new List<Task>();

                    named_double[][] motifs_binary = null;
                    tasks2.Add(Task.Run(() => { motifs_binary = motifs == null ? motifs : motifs.Select(a => a == null ? a : a.Select(b => new named_double() {name = b.name, value = b.value != 0 ? 1 : 0}).ToArray()).ToArray(); }));

                    //named_double[] saac_binary = null;
                    //tasks2.Add(Task.Run(() => { saac_binary = saac == null ? saac : saac.Select(b => new named_double() {name = b.name, value = b.value != 0 ? 1 : 0}).ToArray(); }));

                    named_double[] oaac_binary = null;
                    tasks2.Add(Task.Run(() => { oaac_binary = oaac == null ? oaac : oaac.Select(b => new named_double() {name = b.name, value = b.value != 0 ? 1 : 0}).ToArray(); }));

                    named_double[][] dipeptides_binary = null;
                    tasks2.Add(Task.Run(() => { dipeptides_binary = dipeptides == null ? dipeptides : dipeptides.Select(a => a == null ? a : a.Select(b => new named_double() {name = b.name, value = b.value != 0 ? 1 : 0}).ToArray()).ToArray(); }));

                    Task.WaitAll(tasks2.ToArray<Task>());


                    result.Add((alphabet.id, alphabet.name, motifs, motifs_binary, /*saac, saac_binary,*/ oaac, oaac_binary, average_seq_positions, dipeptides, dipeptides_binary, average_dipeptide_distance));

                    return result;
                });

                a_tasks.Add(a_task);

            }

            Task.WaitAll(a_tasks.ToArray<Task>());

            var ext_result = a_tasks.SelectMany(a => a.Result).ToList();

            return ext_result;
        }


        //public List<(int id, string name, double[] values)> feature_AAOrder(string seq, bool sqrt, bool as_dist)
        //{
        //    var result = new List<(int id, string name, double[] values)>();

        //    var indexes = alphabets.Select(alphabet => alphabet.groups.Select(group => seq.IndexOfAll(group.ToCharArray()).ToList()).ToList()).ToList();

        //    for (var alphabet_index = 0; alphabet_index < alphabets.Count; alphabet_index++)
        //    {
        //        var alphabet = alphabets[alphabet_index];

        //        var alphabet_indexes = indexes[alphabet_index];

        //        var alphabet_pairs = alphabet.groups.SelectMany((a, x) => alphabet.groups.Select((b, y) => (x, y)).ToList()).Where((c, y) => c.x <= c.y).ToList();

        //        var group_result = new double[alphabet_pairs.Count];

        //        for (var alphabet_pairs_index = 0; alphabet_pairs_index < alphabet_pairs.Count; alphabet_pairs_index++)
        //        {
        //            var alphabet_pair = alphabet_pairs[alphabet_pairs_index];
        //            var indexes1 = alphabet_indexes[alphabet_pair.x];
        //            var indexes2 = alphabet_indexes[alphabet_pair.y];

        //            var distances = indexes1.SelectMany(a => indexes2.Select(b => (double)Math.Abs(a - b)).ToList()).ToList();

        //            var value = distances.Average();

        //            value = transform_value(value, seq.Length, sqrt, as_dist);

        //            group_result[alphabet_pairs_index] = value;
        //        }

        //        result.Add((alphabet.id, alphabet.name, group_result));
        //    }

        //    return result;
        //}



        //public List<double[]> feature_AAOrder(string seq, bool sqrt, bool as_dist)
        //{
        //    var result = new List<double[]>();
        //    for (var alphabet_index = 0; alphabet_index < alphabets.Count; alphabet_index++)
        //    {
        //        var alphabet = alphabets[alphabet_index];

        //        var group_result = new double[alphabet.Count];

        //        for (var index = 0; index < alphabet.Count; index++)
        //        {
        //            var member = alphabet[index];

        //            var indexes = member.SelectMany(a => seq.IndexOfAll(a)).ToList();

        //            var distances = indexes.SelectMany(a => indexes.Select(b => (double) Math.Abs(a - b)).ToList()).ToList();

        //            var value = distances.Average();


        //            value = transform_value(value, seq.Length, sqrt, as_dist);

        //            group_result[index] = value;
        //        }

        //        result.Add(group_result);
        //    }

        //    return result;
        //}










        //public static List<(int interval, double value, char aa1, char aa2)> aa_sequence_dipeptide_pairs(int sequence_interval, string seq, aa_alphabet aa_alphabet)
        //{
        //    // get dipeptide counts for set interval

        //    List<string> pairs = null;
        //    if (aa_alphabet == aa_alphabet.Normal) pairs = feature_calcs.aa_pairs_half_with_middle;

        //    var composition = new double[pairs.Count];



        //    for (var index = 0; index < pairs.Count; index++)
        //    {
        //        var pair = pairs[index];
        //        var p1 = pair[0];
        //        var p2 = pair[1];



        //        for (var i = 0; i < seq.Length - sequence_interval; i++)
        //        {
        //            var aa1 = seq[i];
        //            var aa2 = seq[i + sequence_interval];



        //            if ((aa1 == p1 && aa2 == p2) || (aa1 == p2 && aa2 == p1))
        //            {
        //                composition[index]++;
        //            }
        //        }
        //    }

        //    return composition;
        //}


        //public static double aa_sequence_distance(string seq1, string seq2)
        //{
        //    // get distance between 2 sequences by comparing AAs 

        //    bool sqrt = true;
        //    bool as_dist = true;

        //    var distribution1 = aa_distribution(seq1, sqrt, as_dist);
        //    var distribution2 = aa_distribution(seq2, sqrt, as_dist);

        //    var diff = 0d;
        //    for (var i = 0; i < distribution1.Count; i++)
        //    {
        //        diff += Math.Abs(distribution1[i].value - distribution2[i].value);
        //    }

        //    return diff;
        //}


        //public static double aa_group_sequence_distance(string seq1, string seq2)
        //{
        //    // get distance between 2 sequences by comparing AA groups
        //    bool sqrt = true;
        //    bool as_dist = true;

        //    var distribution1 = aa_group_distribution(seq1, sqrt, as_dist);
        //    var distribution2 = aa_group_distribution(seq2, sqrt, as_dist);

        //    var diff = 0d;
        //    for (var i = 0; i < distribution1.Count; i++)
        //    {
        //        diff += Math.Abs(distribution1[i].value - distribution2[i].value); // Math.Abs()?
        //    }

        //    return diff;
        //}



        //public static readonly List<(int index, string group_name, string aa_list)> aa_chem_groups =
        //    new List<(int index, string group_name, string aa_list)>()
        //    {
        //        (0, "Aliphatic", "ILV"),
        //        (1, "Hydroxylic", "TS"),
        //        (2, "Tiny", "AGCS"),
        //        (3, "Small", "VPAGCSTDN"),
        //        (4, "Acidic", "NQ"),
        //        (5, "Positive", "HKR"),
        //        (6, "Negative", "DE"),
        //        (7, "Polar", "CSTNQDEHKRYW"),
        //        (8, "Non-polar", "ILVAMFGP"),
        //        (9, "Charged", "DEHKR"),
        //        (10, "Hydrophobic", "ILVACTMFYWHK"),
        //        (11, "Aromatic", "FYWH"),
        //        (12, "Sulphuric", "MC"),
        //    };
        //public static readonly List<(string group_pair, int index1, int index2/*, List<int> list*/)> aa_chem_group_pairs = aa_chem_groups.SelectMany((a1, x) => aa_chem_groups.Select((a2, y) => (group_pair: $"{a1.group_name}_{a2.group_name}", index1: a1.index, index2: a2.index/*, list: new List<int>()*/)).ToList()).ToList();
        //public static readonly List<(string group_pair, int index1, int index2/*, List<int> list*/)> aa_chem_group_pairs_half_with_middle = aa_chem_groups.SelectMany((a1, x) => aa_chem_groups.Where((a2, y) => x <= y).Select((a2, y) => (group_pair: $"{a1.group_name}_{a2.group_name}", index1: a1.index, index2: a2.index/*, list: new List<int>()*/)).ToList()).ToList();
        ////public static readonly List<(string group_pair, int index1, int index2, List<int> list)> aa_chem_group_pairs_half = aa_chem_groups.SelectMany((a1, x) => aa_chem_groups.Where((a2, y) => x < y).Select((a2, y) => (group_pair: $"{a1.group_name}_{a2.group_name}", index1: a1.index, index2: a2.index, list: new List<int>())).ToList()).ToList();

        //public const string aa_list = "ARNDCQEGHILKMFPSTWYV";

        //public const string group4_list = "1234";
        //public const string group7_list = "1234";

        //public static readonly List<string> aa_pairs = aa_list.SelectMany((a, x) => aa_list.Select((b, y) => $"{a}{b}").ToList()).ToList();
        //public static readonly List<string> aa_pairs_half_with_middle = aa_list.SelectMany((a, x) => aa_list.Where((b, y) => x <= y).Select((b, y) => $"{a}{b}").ToList()).ToList();


        //public static readonly List<string> group4_pairs = group4_list.SelectMany((a, x) => group4_list.Select((b, y) => $"{a}{b}").ToList()).ToList();
        //public static readonly List<string> group4_pairs_half_with_middle = group4_list.SelectMany((a, x) => group4_list.Where((b, y) => x <= y).Select((b, y) => $"{a}{b}").ToList()).ToList();

        //public static readonly List<string> group7_pairs = group7_list.SelectMany((a, x) => group7_list.Select((b, y) => $"{a}{b}").ToList()).ToList();
        //public static readonly List<string> group7_pairs_half_with_middle = group7_list.SelectMany((a, x) => group7_list.Where((b, y) => x <= y).Select((b, y) => $"{a}{b}").ToList()).ToList();

        //public static readonly List<string> aa_pairs_half = aa_list.SelectMany((a, x) => aa_list.Where((b, y) => x < y).Select((b, y) => $"{a}{b}").ToList()).ToList();

        //public static readonly List<string> aa_thrices = aa_list.SelectMany(a => aa_list.SelectMany(b => aa_list.Select(c => ""+a+b+c).ToList()).ToList()).ToList();
        //return pairs as distnace (min,mid,max,average)


        /*
            x y  condition
            a a  x = y
            a b  x < y
            a c  x < y
            b a  x < y  dup
            b b  x = y
            b c  x < y
            c a  x < y  dup
            c b  x < y  dup
            c c  x = y
         */

        //public static List<(char char_value1, char char_value2, List<int> distance_list)> find_distance_all(string text, string chars = aa_list)
        //{
        //    // find distance between all 'chars' within 'text'
        //    var indexes = chars.Select(a => (char_value: a, char_index_list: text.IndexOfAll(a))).ToList();

        //    var distances = indexes.SelectMany((index1, x) => indexes.Where((index2, y) => x <= y).Select((index2, y) => (char_value1: index1.char_value, char_value2: index2.char_value, distance_list: index1.char_index_list.SelectMany(a => index2.char_index_list.Select(b => Math.Abs(b - a))).ToList())).ToList()).ToList();
        //    distances.ForEach(d => d.distance_list.RemoveAll(a => a == 0));

        //    return distances;
        //}

        //public static List<descriptive_stats> aa_group_pair_distribution(string sequence, bool sqrt, bool as_dist)//, bool direction_sensitive = true)
        //{
        //    // ~13*13/2 features

        //    var distances = find_distance_all(sequence);

        //    var pairs1 = feature_calcs.aa_pairs_half_with_middle;

        //    // get list of AA-to-AA distances to convert into PP-to-PP distances
        //    var aa_distances = pairs1.Where(aa_pair => aa_pair[0] <= aa_pair[1]).Select(aa_pair => (aa_pair: aa_pair, distances: distances.FindIndex(b => b.char_value1 == aa_pair[0] && b.char_value2 == aa_pair[1]) == -1 ? new List<int>() : distances.First(b => b.char_value1 == aa_pair[0] && b.char_value2 == aa_pair[1]).distance_list)).ToList();

        //    // chem groups
        //    var pairs2 = feature_calcs.aa_chem_group_pairs_half_with_middle.Select(a => (a.group_pair, a.index1, a.index2, list: new List<int>())).ToList();

        //    foreach (var aa_pair_distance in aa_distances)
        //    {
        //        var groups_where_aa_is_member1 = aa_chem_groups.Where(group => group.aa_list.Contains(aa_pair_distance.aa_pair[0])).ToList();
        //        var groups_where_aa_is_member2 = aa_chem_groups.Where(group => group.aa_list.Contains(aa_pair_distance.aa_pair[1])).ToList();

        //        foreach (var m1 in groups_where_aa_is_member1)
        //        {
        //            foreach (var m2 in groups_where_aa_is_member2)
        //            {
        //                var t = pairs2.Where(a => (a.index1 == m1.index && a.index2 == m2.index) || (a.index1 == m2.index && a.index2 == m1.index)).ToList();

        //                //if (t.Count > 0) // xtodo: is it a xbug to exclude zero values - does it cause differeing number of features (i.e. columns on each row)?
        //                {
        //                    t.ForEach(a => a.list.AddRange(aa_pair_distance.distances));
        //                }
        //            }
        //        }
        //    }

        //    var sv = pairs2.Select(a =>
        //    {
        //        var x = a.list.Select(b => string.IsNullOrWhiteSpace(sequence) ? (double)0 : transform_value(b,sequence.Length,sqrt,as_dist)).ToArray();
        //        return descriptive_stats.get_stat_values(x, a.group_pair);
        //    }).ToList();

        //    return sv;
        //}

        //public static List<descriptive_stats> aa_pair_distribution(string sequence, bool sqrt, bool as_dist)
        //{
        //    // ~20*20/2 features

        //    var distances = find_distance_all(sequence);

        //    var pairs = feature_calcs.aa_pairs_half_with_middle;

        //    //var aa_distances = pairs.Select(aa_pair => (aa_pair:aa_pair, distances: distances.FirstOrDefault(b => (b.char_value1 == aa_pair[0] && b.char_value2 == aa_pair[1])).distance_list)).ToList();
        //    var aa_distances = pairs.Select(aa_pair => (aa_pair: aa_pair, distances: distances.FirstOrDefault(b => (b.char_value1 == aa_pair[0] && b.char_value2 == aa_pair[1]) || (b.char_value1 == aa_pair[1] && b.char_value2 == aa_pair[0])).distance_list)).ToList();

        //    var sv = aa_distances.Select(a =>
        //    {
        //        var x = a.distances.Select(b => string.IsNullOrWhiteSpace(sequence) ? (double)0 : transform_value(b,sequence.Length,sqrt,as_dist)).ToArray();
        //        return descriptive_stats.get_stat_values(x, a.aa_pair);
        //    }).ToList();

        //    return sv;
        //}
    }
}
