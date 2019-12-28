using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading.Tasks;

namespace Coil_DHC_Dataset_Maker
{
    public static class sable
    {
        //public static List<sable_item> sable_data = sable.load(@"c:\betastrands_dataset\sable\OUT_SABLE_RES.txt");


        public class sable_item
        {
            public string seq_name;
            public int seq_index;
            public char amino_acid;
            public char predicted_ss;
            public double predicted_ss_confidence;
            public double prob_h;
            public double prob_e;
            public double prob_c;
            public double relative_burial_value;
            public double relative_burial_confidence;
            public double absolute_burial_value;
            public double entropy_value;
        }

        public static List<sable_item> load(string filename)// = @"c:\betastrands_dataset\sable\OUT_SABLE_RES.txt")
        {
            var lines = File.ReadAllLines(filename).ToList();

            var queries = new List<List<string>>();

            var temp_query = new List<string>();
            for (var i = 0; i <= lines.Count; i++)
            {

                if (i == lines.Count || lines[i].StartsWith("Query:"))
                {
                    if (temp_query != null && temp_query.Count > 0)
                    {
                        queries.Add(temp_query);
                    }

                    temp_query = new List<string>();

                    if (i == lines.Count) break;

                }

                if (!string.IsNullOrWhiteSpace(lines[i]))
                {
                    temp_query.Add(lines[i]);
                }
            }

            var result = new List<sable_item>();

            foreach (var query in queries)
            {
                var query_name = query.First(a => a.StartsWith("Query:")).Split().Last();

                // 3 sectionss
                var SECTION_SS = query.Skip(query.FindIndex(a => a == "SECTION_SS")+1).TakeWhile(a => a != "END_SECTION").SkipWhile(a => !a.StartsWith(">")).Where(a => !a.StartsWith(">")).Select(a => a.Trim()).ToList();

                var SECTION_SS_joined = (
                    query_seq: string.Join("", SECTION_SS.Where((a, i) =>          i % 3 == 0).ToList()).ToCharArray(),
                    predicted_ss_seq: string.Join("", SECTION_SS.Where((a, i) =>   i % 3 == 1).ToList()).ToCharArray(),
                    confidence_level: string.Join("", SECTION_SS.Where((a, i) =>   i % 3 == 2).ToList()).Select(a=>double.Parse(a.ToString()) / 10).ToList()
                    );

                // 4 sections H-> E-> C->
                var SECTION_SS_PROBABILITIES = query.Skip(query.FindIndex(a => a == "SECTION_SS_PROBABILITIES")+1).TakeWhile(a => a != "END_SECTION").SkipWhile(a => !a.StartsWith(">")).Where(a => !a.StartsWith(">")).Select(a => a.Trim()).ToList();

                var SECTION_SS_PROBABILITIES_joined = (
                    query_seq: string.Join("", SECTION_SS_PROBABILITIES.Where((a, i) =>   i % 4 == 0).ToList()).Replace(" ","").ToCharArray(),
                    prob_h: string.Join("", SECTION_SS_PROBABILITIES.Where((a, i) =>   i % 4 == 1).ToList()).Replace("H->", "").Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Select(a => double.Parse(a) / 100).ToList(),
                    prob_e: string.Join("", SECTION_SS_PROBABILITIES.Where((a, i) => i % 4 == 2).ToList()).Replace("E->", "").Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Select(a => double.Parse(a) / 100).ToList(),
                    prob_c: string.Join("", SECTION_SS_PROBABILITIES.Where((a, i) =>  i % 4 == 3).ToList()).Replace("C->", "").Split(new[] { ' ' }, StringSplitOptions.RemoveEmptyEntries).Select(a => double.Parse(a) / 100).ToList()
                    );

                // 3 sections
                var SECTION_SA = query.Skip(query.FindIndex(a => a == "SECTION_SA")+1).TakeWhile(a => a != "END_SECTION").SkipWhile(a => !a.StartsWith(">")).Where(a => !a.StartsWith(">")).Select(a => a.Trim()).ToList();

                var SECTION_SA_joined = (
                    query_seq: string.Join("", SECTION_SA.Where((a, i) =>   i % 3 == 0).ToList()).ToCharArray(),
                    burial_level: string.Join("", SECTION_SA.Where((a, i) =>   i % 3 == 1).ToList()).Select(a => double.Parse(a.ToString()) / 10).ToList(), //0=buried, 9=exposed
                    confidence_level: string.Join("", SECTION_SA.Where((a, i) => i % 3 == 2).ToList()).Select(a => double.Parse(a.ToString()) / 10).ToList()
                );

                //every 2 lines
                var SECTION_SA_ABSOLUTE = query.Skip(query.FindIndex(a => a == "SECTION_SA_ABSOLUTE")+1).TakeWhile(a => a != "END_SECTION").SkipWhile(a => !a.StartsWith(">")).Where(a => !a.StartsWith(">")).Select(a => a.Trim()).ToList();

                var SECTION_SA_ABSOLUTE_joined = (
                    seq: string.Join("", SECTION_SA_ABSOLUTE.Where((a, i) => i % 2 == 0).ToList()).Replace(" ","").ToCharArray(),
                    absolute_burial_value: string.Join(" ", SECTION_SA_ABSOLUTE.Where((a, i) => i % 2 == 1).ToList()).Split(new []{' '},StringSplitOptions.RemoveEmptyEntries).Select(a => double.Parse(a)/100).ToList()
                );

                // every 2 lines "ENTROPY->"
                var SECTION_ENTROPY = query.Skip(query.FindIndex(a => a == "SECTION_ENTROPY")+1).TakeWhile(a => a != "END_SECTION").SkipWhile(a => !a.StartsWith(">")).Where(a => !a.StartsWith(">")).Select(a => a.Trim()).ToList();

                var SECTION_ENTROPY_joined = (
                    seq: string.Join("", SECTION_ENTROPY.Where((a, i) => i % 2 == 0).ToList()).Replace(" ","").ToCharArray(),
                    entropy: string.Join("", SECTION_ENTROPY.Where((a, i) => i % 2 == 1).ToList()).Replace("ENTROPY->","").Split(new char[]{' '},StringSplitOptions.RemoveEmptyEntries).Select(a=>double.Parse(a)).ToList()
                );

                var joined = SECTION_ENTROPY_joined.seq.Select((a, i) => new sable_item()
                {
                    seq_name = query_name,
                    seq_index = i,
                    amino_acid = a,
                    absolute_burial_value = SECTION_SA_ABSOLUTE_joined.absolute_burial_value[i],
                    entropy_value = SECTION_ENTROPY_joined.entropy[i],
                    predicted_ss = SECTION_SS_joined.predicted_ss_seq[i],
                    predicted_ss_confidence = SECTION_SS_joined.confidence_level[i],
                    prob_h = SECTION_SS_PROBABILITIES_joined.prob_h[i],
                    prob_e = SECTION_SS_PROBABILITIES_joined.prob_e[i],
                    prob_c = SECTION_SS_PROBABILITIES_joined.prob_c[i],
                    relative_burial_confidence = SECTION_SA_joined.confidence_level[i],
                    relative_burial_value = SECTION_SA_joined.burial_level[i],

                }).ToList();

                var normalse_ss_prob = true;
                var normalise_entropy = true;
                var normalise_burial_rel = true;
                var normalise_burial_abs = true;

                if (normalse_ss_prob)
                {
                    //var min = joined.Min(a => new double[] { a.prob_h, a.prob_e, a.prob_c }.Min());
                    //var max = joined.Min(a => new double[] { a.prob_h, a.prob_e, a.prob_c }.Max());

                    

                    foreach (var item in joined)
                    {
                        var prob_total = item.prob_h + item.prob_e + item.prob_c;
                        if (prob_total != 0)
                        {
                            item.prob_h = item.prob_h / prob_total;
                            item.prob_e = item.prob_e / prob_total;
                            item.prob_c = item.prob_c / prob_total;
                        }

                       
                    }
                }

                if (normalise_entropy)
                {
                    var min = joined.Min(a => a.entropy_value);
                    var max = joined.Max(a => a.entropy_value);

                }

                if (normalise_burial_rel)
                {
                    var min = joined.Min(a => a.relative_burial_value);
                    var max = joined.Max(a => a.relative_burial_value);

                }

                if (normalise_burial_abs)
                {
                    var min = joined.Min(a => a.absolute_burial_value);
                    var max = joined.Max(a => a.absolute_burial_value);

                }

                result.AddRange(joined);
            }

            return result;

        }
    }
}
