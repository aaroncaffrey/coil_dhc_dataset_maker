using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Coil_DHC_Dataset_Maker
{
    public static class iupred
    {
        public static List<(int index, double short_type_score, double long_type_score, double glob_type_score, double anchor2_score)> load(string pdb_id, char chain_id)
        {
       
            //var seq_index_lookup_file = @"C:\betastrands_dataset\betastrands_dataset_sequences_limited.txt";

            //var seq_lookup_list = File.ReadAllLines(seq_index_lookup_file).ToList();

           // var index = seq_lookup_list.IndexOf(sequence);

            //if (index == -1)

            //{
            //    Console.WriteLine("Warning: iupred2a data missing");
            //    return new List<(int index, double short_type_score, double long_type_score, double glob_type_score, double anchor2_score)>();
            //}

            var iup_files = new List<(string type, string filename)>()
            {
                ("long",$@"C:\betastrands_dataset\iupred2a\{pdb_id}{chain_id}.seq.long.iup"),
                ("short",$@"C:\betastrands_dataset\iupred2a\{pdb_id}{chain_id}.seq.short.iup"),
                ("glob",$@"C:\betastrands_dataset\iupred2a\{pdb_id}{chain_id}.seq.glob.iup"),
            };

            var iup_data = iup_files.SelectMany(a => File.ReadAllLines(a.filename).Where(b=>!b.StartsWith("#") && "0123456789".Contains(b.First())).Select((b,i) =>
            {

                var x= b.Split(new[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);

                var y = (type:a.type, line_index:i, data_index:int.Parse(x[0]), amino_acid:x[1][0], score: double.Parse(x[2]), anchor2:double.Parse(x[3]));

                return y;

            }).ToList()).ToList();

            var iup_data_combined = iup_data.GroupBy(a => a.line_index).Select(a =>
            {
                var list = a.ToList();

                var first = list.First();

                var long_type_score = list.First(b => b.type == "long").score;
                var short_type_score = list.First(b => b.type == "short").score;
                var glob_type_score = list.First(b => b.type == "glob").score;
                var anchor2_score = first.anchor2;

                return (first.line_index, short_type_score, long_type_score, glob_type_score, anchor2_score);

            }).ToList();

            return iup_data_combined;
        }
    }
}
