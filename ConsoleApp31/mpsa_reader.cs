using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace Coil_DHC_Dataset_Maker
{
    public class mpsa_reader
    {
        public string filename = null;
        public string format = null;
        public List<char> ss_column_headers = null;
        public List<mpsa_line_entry> mpsa_matrix = null;

        public List<(char ss, double prob_value, double dist_value)> ss_overall_average = null;
        public List<List<(char ss, double prob_value, double dist_value)>> split_ss_overall_average = null;

        public List<(char ss, int alphabet_id, string alphabet_name, string alphabet_group, double prob_value, double dist_value)> ss_probabilities_per_aa = null;
        public List<List<(char ss, int alphabet_id, string alphabet_name, string alphabet_group, double prob_value, double dist_value)>> split_ss_probabilities_per_aa = null;


        public class mpsa_line_entry
        {
            public mpsa_reader reader = null;
            public int index = -1;
            public char amino_acid = ' ';
            public char predicted_ss_code = ' ';
            public List<char> ss_column_headers = new List<char>() { 'H','E','C','T'};
            public List<(char ss, char amino_acid, double value)> line_prob_values = null;

            public mpsa_line_entry()
            {
                
            }

            public mpsa_line_entry(mpsa_line_entry mpsa_line_entry)
            {
                if (mpsa_line_entry == null) return;

                this.reader = mpsa_line_entry.reader;
                this.index = mpsa_line_entry.index;
                this.amino_acid = mpsa_line_entry.amino_acid;
                this.predicted_ss_code = mpsa_line_entry.predicted_ss_code;
                this.ss_column_headers = mpsa_line_entry.ss_column_headers.ToList();
                this.line_prob_values = mpsa_line_entry.line_prob_values.ToList();
            }

          
        }

        private void mpsa_reader_init()
        {
            var prob_dist = calc_prob_dist(this.mpsa_matrix, this.ss_column_headers);
            this.ss_overall_average = prob_dist.ss_overall_average;
            this.ss_probabilities_per_aa = prob_dist.ss_probabilities_per_aa;

            var mpsa_matrix_split = feature_calcs.split_sequence(mpsa_matrix, 3, 0, false);
            var mpsa_matrix_split_prob_dist = mpsa_matrix_split.Select(a => calc_prob_dist(a, this.ss_column_headers)).ToList();
            this.split_ss_overall_average = mpsa_matrix_split_prob_dist.Select(a => a.ss_overall_average).ToList();
            this.split_ss_probabilities_per_aa = mpsa_matrix_split_prob_dist.Select(a => a.ss_probabilities_per_aa).ToList();


            //var prob_dist = calc_prob_dist(this.mpsa_matrix, this.ss_column_headers);
            //this.ss_overall_average = prob_dist.ss_overall_average;
            //this.ss_probabilities_per_aa = prob_dist.ss_probabilities_per_aa;

            //var mpsa_matrix_split = feature_calcs.split_sequence(mpsa_matrix, 3, true);
            //var mpsa_matrix_split_prob_dist = mpsa_matrix_split.Select(a => calc_prob_dist(a, this.ss_column_headers)).ToList();
            //this.split_ss_overall_average = mpsa_matrix_split_prob_dist.Select(a => a.ss_overall_average).ToList();
            //this.split_ss_probabilities_per_aa = mpsa_matrix_split_prob_dist.Select(a => a.ss_probabilities_per_aa).ToList();


            //var prob_dist = calc_prob_dist(this.mpsa_matrix, this.ss_column_headers);
            //this.ss_overall_average = prob_dist.ss_overall_average;
            //this.ss_probabilities_per_aa = prob_dist.ss_probabilities_per_aa;

            //var mpsa_matrix_split = feature_calcs.split_sequence(mpsa_matrix, 3, true);
            //var mpsa_matrix_split_prob_dist = mpsa_matrix_split.Select(a => calc_prob_dist(a, this.ss_column_headers)).ToList();
            //this.split_ss_overall_average = mpsa_matrix_split_prob_dist.Select(a => a.ss_overall_average).ToList();
            //this.split_ss_probabilities_per_aa = mpsa_matrix_split_prob_dist.Select(a => a.ss_probabilities_per_aa).ToList();


            //var prob_dist = calc_prob_dist(this.mpsa_matrix, this.ss_column_headers);
            //this.ss_overall_average = prob_dist.ss_overall_average;
            //this.ss_probabilities_per_aa = prob_dist.ss_probabilities_per_aa;

            //var mpsa_matrix_split = feature_calcs.split_sequence(mpsa_matrix, 3, true);
            //var mpsa_matrix_split_prob_dist = mpsa_matrix_split.Select(a => calc_prob_dist(a, this.ss_column_headers)).ToList();
            //this.split_ss_overall_average = mpsa_matrix_split_prob_dist.Select(a => a.ss_overall_average).ToList();
            //this.split_ss_probabilities_per_aa = mpsa_matrix_split_prob_dist.Select(a => a.ss_probabilities_per_aa).ToList();
        }

        public mpsa_reader(string format, List<(int index, char amino_acid, char predicted_ss_code, double prob_h, double prob_e, double prob_c)> ss_matrix)
        {
            this.ss_column_headers = new List<char>() {'H', 'E', 'C'};
            this.filename = "";
            this.format = format;
            
            this.mpsa_matrix = ss_matrix.Select(a => new mpsa_line_entry()
            {
                index = a.index,
                amino_acid = a.amino_acid,
                reader = this,
                line_prob_values = new List<(char ss, char amino_acid, double value)>()
                {
                    ('H', a.amino_acid, a.prob_h),
                    ('E',a.amino_acid,a.prob_e),
                    ('C',a.amino_acid,a.prob_c),
                },
                predicted_ss_code = a.predicted_ss_code,
                ss_column_headers = new List<char>() { 'H','E','C' }
            }).ToList();

            mpsa_reader_init();
        }

        public mpsa_reader(mpsa_reader mpsa_reader, List<int> indexes)
        {
            this.filename = mpsa_reader.filename;
            this.format = mpsa_reader.format;
            this.ss_column_headers = mpsa_reader.ss_column_headers.ToList();
            this.mpsa_matrix = mpsa_reader.mpsa_matrix.Where(a=>indexes.Contains(a.index)).Select(a => new mpsa_line_entry(a) { reader = this }).ToList();

            mpsa_reader_init();
        }

        public mpsa_reader(List<mpsa_reader> mpsa_readers)
        {
            this.filename = "consensus";
            this.format = "consensus";
            this.mpsa_matrix = null;
            this.ss_column_headers = "HEC".ToList();
            this.ss_overall_average = null;
            this.ss_probabilities_per_aa = null;

            if (mpsa_readers == null || mpsa_readers.Count == 0) return;

            var matrixes = mpsa_readers.Where(a=>a.mpsa_matrix!=null&&a.mpsa_matrix.Count>0).SelectMany(a => a.mpsa_matrix).ToList();

            this.mpsa_matrix = matrixes.GroupBy(a => a.index).Select(a =>
            {
                var index = a.Key;
                var lines = a.ToList();

                var line_prob_values = lines
                    .SelectMany(b => b.line_prob_values.Where(c=>"HEC".Contains(c.ss)))
                    .GroupBy(b => b.ss)
                    .Select(b => (ss: b.Key, amino_acid: b.First().amino_acid, value: b.Select(c => c.value).Average())).ToList();

                var predicted_ss_code = line_prob_values.OrderByDescending(b => b.value).First().ss;

                var x = new mpsa_line_entry()
                {
                    amino_acid = lines[0].amino_acid,
                    reader = this,
                    index = index,
                    predicted_ss_code = predicted_ss_code,
                    line_prob_values = line_prob_values,
                    ss_column_headers = "HEC".ToList()
                };
                return x;
            }).ToList();


            mpsa_reader_init();
        }

        //public static List<(char predicted_ss, char amino_acid, List<(char ss, char amino_acid, double average_value)> values)> consensus_individual_probability(List<mpsa_reader> mpsa_readers, List<int> indexes)
        //{
        //    mpsa_readers = mpsa_readers.Where(a => a.mpsa_matrix != null && a.mpsa_matrix.Count > 0).ToList();

        //    mpsa_readers = mpsa_readers.Select(a => new mpsa_reader(a, indexes)).ToList();

        //    var x0 = mpsa_readers.SelectMany(a => a.mpsa_matrix.Where(b=>"HECT".Contains(b.predicted_ss_code)).ToList()).GroupBy(a => a.index).ToList();
        //    var x1 = x0.Select(a => a.SelectMany(b => b.line_prob_values.Where(c=>"HECT".Contains(c.ss)).ToList()).ToList()).ToList();
        //    var x2 = x1.Select(a => a.GroupBy(b => b.ss).Select(b => (ss:b.Key, b.First().amino_acid, average_value: b.Select(c => c.value).DefaultIfEmpty(0).Average())).ToList()).ToList();

        //    var x3 = x2.Select(a => (a.OrderByDescending(b=>b.average_value).First().ss,a.First().amino_acid,a)).ToList();

        //    return x3;
        //}

        //public static List<(char ss, double average_value)> consensus_overall_probability(List<mpsa_reader> mpsa_readers, List<int> indexes)
        //{
        //    mpsa_readers = mpsa_readers.Where(a => a.mpsa_matrix != null && a.mpsa_matrix.Count > 0).ToList();

        //    mpsa_readers = mpsa_readers.Select(a => new mpsa_reader(a, indexes)).ToList();

        //    var calculate_from_individual = true;

        //    if (calculate_from_individual)
        //    {
        //        var y = consensus_individual_probability(mpsa_readers, indexes);

        //        var x = y.SelectMany(a => a.values.Where(b => "HECT".Contains(b.ss)).ToList()).GroupBy(a => a.ss).Select(a => (ss: a.Key, averge_value: a.Select(b => b.average_value).DefaultIfEmpty(0).Average())).OrderByDescending(a => a.averge_value).ToList();

        //        return x;
        //    }
        //    else
        //    {

        //        var x = mpsa_readers.SelectMany(a => a.ss_overall_average).GroupBy(a => a.ss).Select(a => (ss: a.Key, average: a.Select(b => b.prob_value).DefaultIfEmpty(0).Average())).OrderByDescending(a => a.average).ToList();

        //        return x;
        //    }
        //}

        //public static List<(char ss, double average_value)> consensus_overall_distribution
        //    (List<mpsa_reader> mpsa_readers, List<int> indexes)
        //{
        //    mpsa_readers = mpsa_readers.Where(a => a.mpsa_matrix != null && a.mpsa_matrix.Count > 0).ToList();

        //    mpsa_readers = mpsa_readers.Select(a => new mpsa_reader(a, indexes)).ToList();

        //    var as_dist = true;
        //    var sqrt = false;

           
        //        var y = consensus_individual_probability(mpsa_readers, indexes);

        //        var ss_seq = y.Select(a => a.predicted_ss).ToList();

        //        var x = "HECT".Select(ss => (ss,feature_calcs.transform_value(ss_seq.Count(b => ss == b),ss_seq.Count,sqrt,as_dist))).ToList();
        //        return x;
            
        //}

        private List<mpsa_line_entry> mpsa_matrix_fake(string format, List<char> matrix_column_headers, List<char> amino_acids = null)
        {
            var mpsa_matrix = new List<mpsa_line_entry>();

            for (var index = 0; index < amino_acids.Count; index++)
            {
                mpsa_matrix.Add(new mpsa_line_entry()
                {
                    ss_column_headers = matrix_column_headers,
                    reader = this,
                    amino_acid = amino_acids[index],
                    index = index,
                    line_prob_values = matrix_column_headers.Select(a => (a, amino_acids[index], 0d)).ToList(),
                    predicted_ss_code = ' '
                });
            }


            return mpsa_matrix;
        }

        public mpsa_reader(string filename, List<char> amino_acids = null)
        {
            this.filename = filename;
            this.format = Path.GetExtension(filename).Substring(1).ToLowerInvariant();

            this.mpsa_matrix = null;
            this.ss_overall_average = null;
            this.ss_probabilities_per_aa = null;

            if (!File.Exists(this.filename) || new FileInfo(this.filename).Length == 0)
            {
                var matrix_column_headers = mpsa_reader.secondary_structure_codes.First(a => a.format == this.format).ss_codes.ToList();
                this.ss_column_headers = matrix_column_headers;
                this.mpsa_matrix = mpsa_matrix_fake(format, matrix_column_headers, amino_acids);
            }
            else
            {
                var data = read_mpsa_file(filename);

                for (var i = 0; i < data.individual_residue_predictions.Count; i++)
                {
                    data.individual_residue_predictions[i].reader = this;
                }

                this.ss_column_headers = data.matrix_colum_headers;
                this.mpsa_matrix = data.individual_residue_predictions;

                if (this.mpsa_matrix == null || this.mpsa_matrix.Count == 0)
                {
                    //throw new Exception();
                    this.mpsa_matrix = mpsa_matrix_fake(format, data.matrix_colum_headers, amino_acids);
                }

                if (amino_acids != null && mpsa_matrix != null && this.mpsa_matrix.Count != amino_acids.Count) throw new Exception();

                mpsa_reader_init();
            }
        }

        public static (List<(char ss, double, double)> ss_overall_average, List<(char ss, int alphabet_id, string alphabet_name, string alphabet_group, double prob_value, double dist_value)> ss_probabilities_per_aa) calc_prob_dist(List<mpsa_line_entry> mpsa_matrix, List<char> ss_column_headers)
        {
            var as_dist = true;
            var sqrt = false;

            var ss_codes = "HEC";//"HECT";

            var seq_length = (double)mpsa_matrix.Count(b => ss_codes.Contains(b.predicted_ss_code));

            var all_probability_data  // each AA in sequence has 3 or 4 lines
                = ss_column_headers.Where(a => ss_codes.Contains(a)).SelectMany(ss =>
                    (mpsa_matrix == null || mpsa_matrix.Count == 0) ? new List<(char ss, char amino_acid, double value)>() :
                        mpsa_matrix.Select(b => (ss, b.amino_acid, b.line_prob_values.First(c => c.ss == ss).value)).ToList()).ToList();

            // 4 SS
            var all_distribution_average = ss_column_headers.Where(a => ss_codes.Contains(a))
                .Select(ss => (mpsa_matrix == null || mpsa_matrix.Count == 0) ?
                (ss:ss, average:0d) :
                (ss:ss, average:feature_calcs.transform_value((double)mpsa_matrix.Count(a => a.predicted_ss_code == ss), seq_length, sqrt, as_dist))).ToList();

            var all_probability_data_grouped = all_probability_data.GroupBy(a => a.ss).ToList();

            // 4 SS average probability
            // note: probability should not be divded by sequence length because averaging does that anyway
            var all_probabilities_average = all_probability_data_grouped.Select(a => (ss:a.Key, average:a.Select(b => b.value).DefaultIfEmpty(0).Average())).ToList();


            List<(char ss, double, double)> ss_overall_average = ss_column_headers.Where(a => ss_codes.Contains(a)).Select(ss => (ss, all_probabilities_average.FirstOrDefault(a => a.ss==ss).average, all_distribution_average.FirstOrDefault(a => a.ss==ss).average)).ToList();

            // 4 SS x 20 AA+ all other alphabets
            //var all_distribution_per_aa = new List<(char ss, int alphabet_id, string alphabet_group, double value)>();

            //var all_probabilities_per_aa = new List<(char ss, int alphabet_id, string alphabet_group, double value)>();

            List<(char ss, int alphabet_id, string alphabet_name, string alphabet_group, double prob_value, double dist_value)> ss_probabilities_per_aa = new List<(char ss, int alphabet_id, string alphabet_name, string alphabet_group, double prob_value, double dist_value)>();
            foreach (var alphabet in feature_calcs.aa_alphabets)
            {
                foreach (var group in alphabet.groups)
                {
                    var dist = ss_column_headers.Where(a => ss_codes.Contains(a)).Select(ss => (mpsa_matrix == null || mpsa_matrix.Count == 0) ?
                        (ss, alphabet.id, group, average: 0d) :
                        (ss, alphabet.id, group, average: feature_calcs.transform_value(
                            value: (double)mpsa_matrix.Count(a => group.group_amino_acids.Contains(a.amino_acid) && a.predicted_ss_code == ss),
                            length: (double)seq_length, sqrt, as_dist))).ToList();

                    //this.all_distribution_per_aa.AddRange(x);
            
                    var prob = all_probability_data_grouped.Select(a => (ss:a.Key, alphabet.id, alphabet.name, group, average: a.Where(b => group.group_amino_acids.Contains(b.amino_acid)).Select(b => b.value).DefaultIfEmpty(0).Average())).ToList();
                   // this.all_probabilities_per_aa.AddRange(xx);
                    ss_probabilities_per_aa.AddRange(ss_column_headers.Where(a => ss_codes.Contains(a)).Select(ss=>(ss,alphabet.id,alphabet.name, group.group_amino_acids, prob.FirstOrDefault(a => a.ss == ss).average, dist.FirstOrDefault(a => a.ss == ss).average)).ToList());
                    //Console.WriteLine();
                }
            }

            return (ss_overall_average, ss_probabilities_per_aa);
            //Console.WriteLine();
        }


        public static List<(int id, string format, string ss_codes)> secondary_structure_codes = new List<(int id, string format, string ss_codes)>()
        {
            (00, "dpm", "HETCF"), // extra: T F
            (01, "dsc", "HEC"), // non-standard format
            (02, "gor1", "HETC"), // extra: T
            (03, "gor3", "HEC"),
            (04, "hnn", "HEC"),
            (05, "mpsa", "HEC"),
            (06, "phd", "HEC"/*"HECX"*/), // extra: X (reliability) // non-standard format
            (07, "preda", "HEC"), // non-standard format
            (08, "sopm", "HETC"), // extra: T
            (09, "psipred-swissprot", "CHE"),
            (10, "psipred-uniref90", "CHE"),
        };

        //HEC

        public static void mpsa_normalise_matrix_all(List<mpsa_reader> matrices)
        {
            var all_values = matrices.SelectMany(a => a.mpsa_matrix.SelectMany(b => b.line_prob_values.Select(c => c.value).ToList()).ToList()).ToList();
            var all_min = all_values.Min();
            var all_max = all_values.Max();

            foreach (var matrix in matrices)
            {
                for (var i = 0; i < matrix.mpsa_matrix.Count; i++)
                {
                    for (var j = 0; j < matrix.mpsa_matrix[i].line_prob_values.Count; j++)
                    {
                        if ((all_max - all_min) == 0) matrix.mpsa_matrix[i].line_prob_values[j] = (matrix.mpsa_matrix[i].line_prob_values[j].ss, matrix.mpsa_matrix[i].line_prob_values[j].amino_acid, 0);
                        else matrix.mpsa_matrix[i].line_prob_values[j] = (matrix.mpsa_matrix[i].line_prob_values[j].ss, matrix.mpsa_matrix[i].line_prob_values[j].amino_acid, (matrix.mpsa_matrix[i].line_prob_values[j].value - all_min) / (all_max - all_min));
                    }
                }
            }
        }

        public void mpsa_normalise_matrix()
        {
            var all_values = mpsa_matrix.SelectMany(a => a.line_prob_values.Select(b => b.value).ToList()).Distinct().ToList();
            var all_min = all_values.Min();
            var all_max = all_values.Max();

            for (var i = 0; i < mpsa_matrix.Count; i++)
            {
                for (var j = 0; j < mpsa_matrix[i].line_prob_values.Count; j++)
                {
                    if ((all_max - all_min) == 0) mpsa_matrix[i].line_prob_values[j] = (mpsa_matrix[i].line_prob_values[j].ss, mpsa_matrix[i].line_prob_values[j].amino_acid, 0);
                    else mpsa_matrix[i].line_prob_values[j] = (mpsa_matrix[i].line_prob_values[j].ss, mpsa_matrix[i].line_prob_values[j].amino_acid, (mpsa_matrix[i].line_prob_values[j].value - all_min) / (all_max - all_min));
                }
            }
        }

        public void mpsa_normalise_rows()
        {
            for (var i = 0; i < mpsa_matrix.Count; i++)
            {
                var row_values = mpsa_matrix[i].line_prob_values.Select(a => a.value).ToList();
                var row_min = row_values.Min();
                var row_max = row_values.Max();

                for (var j = 0; j < mpsa_matrix[i].line_prob_values.Count; j++)
                {
                    if ((row_max - row_min) == 0) mpsa_matrix[i].line_prob_values[j] = (mpsa_matrix[i].line_prob_values[j].ss, mpsa_matrix[i].line_prob_values[j].amino_acid, 0);
                    else mpsa_matrix[i].line_prob_values[j] = (mpsa_matrix[i].line_prob_values[j].ss, mpsa_matrix[i].line_prob_values[j].amino_acid, (mpsa_matrix[i].line_prob_values[j].value - row_min) / (row_max - row_min));
                }
            }
        }

        public void mpsa_normalise_columns()
        {
            for (var h = 0; h < ss_column_headers.Count; h++)
            {
                var col_values = mpsa_matrix.Select(a => a.line_prob_values.First(b => b.ss == ss_column_headers[h]).value).ToList();
                var col_min = col_values.Min();
                var col_max = col_values.Max();

                for (var i = 0; i < mpsa_matrix.Count; i++)
                {
                    var col = mpsa_matrix[i].line_prob_values.FindIndex(a => a.ss == ss_column_headers[h]);

                    if ((col_max - col_min) == 0) mpsa_matrix[i].line_prob_values[col] = (mpsa_matrix[i].line_prob_values[col].ss, mpsa_matrix[i].line_prob_values[col].amino_acid, 0);
                    else mpsa_matrix[i].line_prob_values[col] = (mpsa_matrix[i].line_prob_values[col].ss, mpsa_matrix[i].line_prob_values[col].amino_acid, (mpsa_matrix[i].line_prob_values[col].value - col_min) / (col_max - col_min));
                }
            }
        }

        private static (List<char> matrix_colum_headers, List<mpsa_line_entry> individual_residue_predictions) read_mpsa_file(string filename)
        {
            // common MPSA matrix format:
            // DPM, GOR1, GOR2, MPSA, HNN, SOPM

            var convert_t_to_c = true;
            var ignore_f = true;
            var ignore_x = true;

            var format = Path.GetExtension(filename).Substring(1);

            var mpsa_extensions_supported = new List<string> { "dpm", "gor1", "gor2", "gor3", "mpsa", "hnn", "sopm", "psipred-uniref90", "psipred-swissprot" };
            var other_extensions_supported = new List<string> { "dsc", "phd", "preda" };
            mpsa_extensions_supported.AddRange(other_extensions_supported);

            if (!mpsa_extensions_supported.Contains(format, StringComparer.InvariantCultureIgnoreCase))
            {
                Console.WriteLine("File format not supported: " + filename);
                return (null, null);
            }

            var matrix_colum_headers = mpsa_reader.secondary_structure_codes.First(a => a.format == format).ss_codes.ToList();

            var lines = new List<string>();

            if (File.Exists(filename) && new FileInfo(filename).Length > 0)
            {
                lines = File.ReadAllLines(filename).ToList();
            }
            

            if (format.Equals("preda",StringComparison.InvariantCultureIgnoreCase))
            {
                var start_marker = lines.FindIndex(a => a.StartsWith("HEADER  |- Residue -|  Pred  Rel      NAli   Asn"));
                lines = lines.Skip(start_marker + 1).ToList();
                lines = lines.Where(a => a.StartsWith("PRED")).ToList();

                var check_matrix_colum_headers = "HEC".ToList();

                if (!matrix_colum_headers.SequenceEqual(check_matrix_colum_headers)) throw new Exception();

                var parsed_matrix1 = new List<mpsa_line_entry>();

                //foreach (var line in lines)
                for (var index = 0; index < lines.Count; index++)
                {
                    var s = lines[index].Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                    //0:PRED  1:216    2:GLY    3:G  4:c     5:0.000    6:0      7:?

                    var values = new List<(char ss, char amino_acid, double value)>();
                    //var index = int.Parse(s[1]);
                    var amino_acid = s[3].ToUpperInvariant()[0];
                    var predicted_ss_code = s[4].ToUpperInvariant()[0];
                    var reliability = double.Parse(s[5]);

                    if (reliability < 0.445d) reliability = 0.445d;

                    var rel_remaining = (1d - reliability) / 2d;

                    if (!matrix_colum_headers.Contains(predicted_ss_code)) throw new Exception();

                    foreach (var h in matrix_colum_headers)
                    {
                        if (h == predicted_ss_code)
                        {
                            values.Add((h, amino_acid, reliability));
                        }
                        else
                        {
                            values.Add((h, amino_acid, rel_remaining));
                        }
                    }

                    var h_exists = matrix_colum_headers.Contains('H');
                    var e_exists = matrix_colum_headers.Contains('E');
                    var c_exists = matrix_colum_headers.Contains('C');
                    var t_exists = matrix_colum_headers.Contains('T');
                    var f_exists = matrix_colum_headers.Contains('F');

                    var prob_h = !h_exists ? 0d : values[values.FindIndex(a => a.ss == 'H')].value;
                    var prob_e = !e_exists ? 0d : values[values.FindIndex(a => a.ss == 'E')].value;
                    var prob_c = !c_exists ? 0d : values[values.FindIndex(a => a.ss == 'C')].value;
                    var prob_t = !t_exists ? 0d : values[values.FindIndex(a => a.ss == 'T')].value;
                    var prob_f = !f_exists ? 0d : values[values.FindIndex(a => a.ss == 'F')].value;

                    if (convert_t_to_c && t_exists)
                    {
                        t_exists = false;
                        prob_c += prob_t;
                        prob_t = 0d;
                    }

                    if (ignore_f && f_exists)
                    {
                        f_exists = false;
                        prob_f = 0d;
                    }

                    var prob_adjust_values = new List<double>();
                    if (h_exists) prob_adjust_values.Add(prob_h);
                    if (e_exists) prob_adjust_values.Add(prob_e);
                    if (c_exists) prob_adjust_values.Add(prob_c);
                    if (t_exists) prob_adjust_values.Add(prob_t);
                    //if (f_exists) prob_adjust_values.Add(prob_f);

                    var prob_adjust_range = prob_adjust_values.Min();

                    if (prob_adjust_range < 0)
                    {
                        var prob_adjust_range_abs = Math.Abs(prob_adjust_range);
                        if (h_exists) prob_h = prob_h + prob_adjust_range_abs;
                        if (e_exists) prob_e = prob_e + prob_adjust_range_abs;
                        if (c_exists) prob_c = prob_c + prob_adjust_range_abs;
                        if (t_exists) prob_t = prob_t + prob_adjust_range_abs;
                        //if (f_exists) prob_f = prob_f + prob_adjust_range_abs;
                    }

                    var prob_total = 0.0;
                    if (h_exists) prob_total += prob_h;
                    if (e_exists) prob_total += prob_e;
                    if (c_exists) prob_total += prob_c;
                    if (t_exists) prob_total += prob_t;
                   // if (f_exists) prob_total += prob_f;

                    if (prob_total != 0)
                    {
                        if (h_exists) prob_h = prob_h / prob_total;
                        if (e_exists) prob_e = prob_e / prob_total;
                        if (c_exists) prob_c = prob_c / prob_total;
                        if (t_exists) prob_t = prob_t / prob_total;
                        //if (f_exists) prob_f = prob_f / prob_total;
                    }

                    values = new List<(char ss, char amino_acid, double value)>();
                    if (h_exists) values.Add(('H', amino_acid, prob_h));
                    if (e_exists) values.Add(('E', amino_acid, prob_e));
                    if (c_exists) values.Add(('C', amino_acid, prob_c));
                    if (t_exists) values.Add(('T', amino_acid, prob_t));
                    //if (f_exists) values.Add(('F', amino_acid, prob_f));
                    
                   

                    var mpsa_entry = new mpsa_line_entry()
                    {
                        
                        index = index,
                        amino_acid = amino_acid,
                        ss_column_headers = matrix_colum_headers.Where(a => a != 'F' && (!convert_t_to_c || a != 'T')).ToList(),
                        predicted_ss_code = predicted_ss_code,
                        line_prob_values = values
                    };

                    parsed_matrix1.Add(mpsa_entry);

                }

                return (matrix_colum_headers, parsed_matrix1);
            }

            if (format.Equals("dsc",StringComparison.InvariantCultureIgnoreCase))
            {
                var start_marker = lines.FindIndex(a => a.StartsWith("NO.  RES   DSC_SEC PROB_H    PROB_E    PROB_C"));
                lines = lines.Skip(start_marker + 1).ToList();
                lines = lines.Where(a => !string.IsNullOrWhiteSpace(a)).ToList();

                //NO.  RES   DSC_SEC PROB_H    PROB_E    PROB_C

                var check_matrix_colum_headers = "HEC".ToList();

                if (!matrix_colum_headers.SequenceEqual(check_matrix_colum_headers)) throw new Exception();


                var parsed_matrix1 = new List<mpsa_line_entry>();

                for (var index = 0; index < lines.Count; index++)
                {
                    var s = lines[index].Split(new char[] { ' ', '\t' }, StringSplitOptions.RemoveEmptyEntries);
                    //NO.  RES   DSC_SEC PROB_H    PROB_E    PROB_C
                    //0:1     1:a      2:C     3:0.070     4:0.040     5:0.890     

                    //var index = int.Parse(s[0]);
                    var amino_acid = s[1].ToUpperInvariant()[0];
                    var predicted_ss_code = s[2][0];
                    var prob_h = double.Parse(s[3]);
                    var prob_e = double.Parse(s[4]);
                    var prob_c = double.Parse(s[5]);

                  
                    var prob_adjust_range = new double[] { prob_h, prob_e, prob_c }.Min();

                    if (prob_adjust_range < 0)
                    {
                        var prob_adjust_range_abs = Math.Abs(prob_adjust_range);
                        prob_h = prob_h + prob_adjust_range_abs;
                        prob_e = prob_e + prob_adjust_range_abs;
                        prob_c = prob_c + prob_adjust_range_abs;
                    }

                    var prob_total = prob_h + prob_e + prob_c;
                    if (prob_total != 0)
                    {
                        prob_h = prob_h / prob_total;
                        prob_e = prob_e / prob_total;
                        prob_c = prob_c / prob_total;
                    }

                    var values = new List<(char ss, char amino_acid, double value)>() { ('H', amino_acid, prob_h), ('E', amino_acid, prob_e), ('C', amino_acid, prob_c), };

                    if (!matrix_colum_headers.Contains(predicted_ss_code)) throw new Exception();

                    var mpsa_entry = new mpsa_line_entry()
                    {
                        index = index,
                        amino_acid = amino_acid,
                        ss_column_headers = matrix_colum_headers.Where(a => a != 'F' && (!convert_t_to_c || a != 'T')).ToList(),
                        predicted_ss_code = predicted_ss_code,
                        line_prob_values = values
                    };

                    parsed_matrix1.Add(mpsa_entry);
                }


                return (matrix_colum_headers, parsed_matrix1);

            }

            if (format.Equals("phd",StringComparison.InvariantCultureIgnoreCase))
            {
                lines = lines.Where(a => !string.IsNullOrWhiteSpace(a) && a.Length > 9 && !a.StartsWith("*") && !a.StartsWith("END")).ToList();
                lines = lines.Select(a => a.Substring(9)).ToList();

                var titles = lines.Select(a => a.Substring(0, 4).Replace("-", "").Trim().ToUpperInvariant()).Distinct().Where(a => !string.IsNullOrWhiteSpace(a)).ToList();
                var merged_lines = new List<(string title, string line, int len)>();
                foreach (var t in titles)
                {
                    merged_lines.Add((t, "", 0));
                }

                foreach (var line in lines)
                {
                    var line_title = line.Substring(0, 4).Replace("-", "").Trim().ToUpperInvariant();
                    if (string.IsNullOrWhiteSpace(line_title)) continue;

                    var line_data = line.Substring(5).Replace("|", "");

                    if (line_title == "PHD")
                    {
                        line_data = line_data.Replace(" ", "C").Replace("L", "C");
                    }

                    var title_index = merged_lines.FindIndex(a => a.title == line_title);
                    var merged_line = merged_lines[title_index].line + line_data;

                    merged_lines[title_index] = (line_title, merged_line, merged_line.Length);
                }

                if (merged_lines.Select(a => a.line.Length).Distinct().Count() != 1) throw new Exception();

                var aa_line = merged_lines.First(a => a.title == "AA").line;
                var rel_line = merged_lines.First(a => a.title == "REL").line;
                var phd_line = merged_lines.First(a => a.title == "PHD").line;

                var prob_h_line = merged_lines.First(a => a.title == "PRH").line;
                var prob_e_line = merged_lines.First(a => a.title == "PRE").line;
                var prob_l_line = merged_lines.First(a => a.title == "PRL").line;


                var check_matrix_colum_headers = "HEC".ToList();//"HECX".ToList();

                if (!matrix_colum_headers.SequenceEqual(check_matrix_colum_headers)) throw new Exception();

                var parsed_matrix1 = new List<mpsa_line_entry>();

                for (var index = 0; index < aa_line.Length; index++)
                {
                    //var seq_index = index + 1;
                    var amino_acid = aa_line[index];
                    var predicted_ss_code = phd_line[index];
                    if (predicted_ss_code == 'L' || predicted_ss_code == ' ')
                    {
                        predicted_ss_code = 'C';
                    }

                    var reliability = double.Parse("" + rel_line[index]);
                    var prob_h = double.Parse("" + prob_h_line[index]);
                    var prob_e = double.Parse("" + prob_e_line[index]);
                    var prob_c = double.Parse("" + prob_l_line[index]);

                    var prob_adjust_range = new double[] {prob_h, prob_e, prob_c}.Min();

                    if (prob_adjust_range < 0)
                    {
                        var prob_adjust_range_abs = Math.Abs(prob_adjust_range);
                        prob_h = prob_h + prob_adjust_range_abs;
                        prob_e = prob_e + prob_adjust_range_abs;
                        prob_c = prob_c + prob_adjust_range_abs;
                    }

                    var prob_total = prob_h + prob_e + prob_c;
                    if (prob_total != 0)
                    {
                        prob_h = prob_h / prob_total;
                        prob_e = prob_e / prob_total;
                        prob_c = prob_c / prob_total;
                    }

                    var values = new List<(char ss, char amino_acid, double value)>();
                    values.Add(('H', amino_acid, prob_h));
                    values.Add(('E', amino_acid, prob_e));
                    values.Add(('C', amino_acid, prob_c));
                    if (!ignore_x) values.Add(('X', amino_acid, reliability));

                    if (!matrix_colum_headers.Contains(predicted_ss_code)) throw new Exception();


                    var mpsa_entry = new mpsa_line_entry()
                    {
                        index = index,
                        amino_acid = amino_acid,
                        ss_column_headers = matrix_colum_headers.Where(a => (!ignore_f || a != 'F') && (!convert_t_to_c || a != 'T') && (!ignore_x || a != 'X')).ToList(),
                        predicted_ss_code = predicted_ss_code,
                        line_prob_values = values
                    };

                    parsed_matrix1.Add(mpsa_entry);

                }

                return (matrix_colum_headers, parsed_matrix1);
            }

            var psipred_marker = lines.FindIndex(a=> a.StartsWith("# PSIPRED VFORMAT"));//# PSIPRED VFORMAT (PSIPRED V4.0)

            var mpsa_code_marker = lines.FindIndex(a => a.StartsWith("MPSA code :"));

            var matrix_start = mpsa_code_marker + 3; // +2 or +3 ?

            if (psipred_marker > -1) matrix_start = psipred_marker + 2;

            
            var matrix_rows = 0;
            var matrix_columns = 0;



            if (mpsa_code_marker > -1)
            {
                var matrix_size_and_fields = lines[matrix_start - 1].Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries);
                matrix_rows = int.Parse(matrix_size_and_fields.First());
                matrix_columns = int.Parse(matrix_size_and_fields.Skip(1).First());
                var check_matrix_colum_headers1 = matrix_size_and_fields.Skip(2).First().Select(a => a).ToList();
                if (!matrix_colum_headers.SequenceEqual(check_matrix_colum_headers1)) throw new Exception();
            }
            else if (psipred_marker > -1)
            {
                matrix_rows = lines.Count - 2;
                matrix_columns = 3;
            }



            var matrix = lines.Skip(matrix_start).Take(matrix_rows).ToList();

            var parsed_matrix = new List<mpsa_line_entry>();
            for (var index = 0; index < matrix.Count; index++)
            {
                var line = matrix[index];
                var s = line.Split(new char[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).ToList();

                var predicted_ss = ' ';
                if (mpsa_code_marker> -1 ) predicted_ss = s.First()[0];
                else if (psipred_marker > -1) predicted_ss = s[2][0];

                if (!matrix_colum_headers.Contains(predicted_ss)) throw new Exception();

                var amino_acid = ' ';
                if (mpsa_code_marker>-1) amino_acid=s.Last()[0];
                else if (psipred_marker > -1) amino_acid = s[1][0];

                if (mpsa_code_marker>-1) s = s.Skip(1).Take(matrix_columns).ToList();
                else if (psipred_marker > -1) s = s.Skip(3).Take(matrix_columns).ToList();

                var values = s.Select((a, i) => (ss: matrix_colum_headers[i], amino_acid:amino_acid, value: double.Parse(a))).ToList();
                values = values.Where(a => a.ss != 'F').ToList();
                
                var h_exists = matrix_colum_headers.Contains('H');
                var e_exists = matrix_colum_headers.Contains('E');
                var c_exists = matrix_colum_headers.Contains('C');
                var t_exists = matrix_colum_headers.Contains('T');
                //var f_exists = matrix_colum_headers.Contains('F');

                var prob_h = !h_exists ? 0d : values[values.FindIndex(a => a.ss == 'H')].value;
                var prob_e = !e_exists ? 0d : values[values.FindIndex(a => a.ss == 'E')].value;
                var prob_c = !c_exists ? 0d : values[values.FindIndex(a => a.ss == 'C')].value;
                var prob_t = !t_exists ? 0d : values[values.FindIndex(a => a.ss == 'T')].value;
                //var prob_f = !f_exists ? 0d : values[values.FindIndex(a => a.ss == 'F')].value;

                if (convert_t_to_c && t_exists)
                {
                    t_exists = false;
                    prob_c += prob_t;
                    prob_t = 0d;
                }

                //if (ignore_f && f_exists)
                //{
                //    f_exists = false;
                //    prob_f = 0d;
                //}

                var prob_adjust_values = new List<double>();
                if (h_exists) prob_adjust_values.Add(prob_h);
                if (e_exists) prob_adjust_values.Add(prob_e);
                if (c_exists) prob_adjust_values.Add(prob_c);
                if (t_exists) prob_adjust_values.Add(prob_t);
                //if (f_exists) prob_adjust_values.Add(prob_f);

                var prob_adjust_range = prob_adjust_values.Min();

                if (prob_adjust_range < 0)
                {
                    var prob_adjust_range_abs = Math.Abs(prob_adjust_range);
                    if (h_exists) prob_h = prob_h + prob_adjust_range_abs;
                    if (e_exists) prob_e = prob_e + prob_adjust_range_abs;
                    if (c_exists) prob_c = prob_c + prob_adjust_range_abs;
                    if (t_exists) prob_t = prob_t + prob_adjust_range_abs;
                    //if (f_exists) prob_f = prob_f + prob_adjust_range_abs;
                }

                var prob_total = 0.0;
                if (h_exists) prob_total += prob_h;
                if (e_exists) prob_total += prob_e;
                if (c_exists) prob_total += prob_c;
                if (t_exists) prob_total += prob_t;
                //if (f_exists) prob_total += prob_f;

                if (prob_total != 0)
                {
                    if (h_exists) prob_h = prob_h / prob_total;
                    if (e_exists) prob_e = prob_e / prob_total;
                    if (c_exists) prob_c = prob_c / prob_total;
                    if (t_exists) prob_t = prob_t / prob_total;
                    //if (f_exists) prob_f = prob_f / prob_total;
                }

                values = new List<(char ss, char amino_acid, double value)>();
                if (h_exists) values.Add(('H', amino_acid, prob_h));
                if (e_exists) values.Add(('E', amino_acid, prob_e));
                if (c_exists) values.Add(('C', amino_acid, prob_c));
                if (t_exists) values.Add(('T', amino_acid, prob_t));
                //if (f_exists) values.Add(('F', amino_acid, prob_f));
                
                var mpsa_entry = new mpsa_line_entry()
                {
                    index = index,
                    amino_acid = amino_acid,
                    ss_column_headers = matrix_colum_headers.Where(a=>a!='F' && (!convert_t_to_c || a!='T')).ToList(),
                    predicted_ss_code = predicted_ss,
                    line_prob_values = values
                };

                parsed_matrix.Add(mpsa_entry);
            }

            return (matrix_colum_headers, parsed_matrix);
        }
    }
}
