using System;
using System.CodeDom;
using System.Collections.Generic;
using System.Diagnostics;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Net;
using System.Runtime.InteropServices;
using System.Threading;
using System.Threading.Tasks;
using RDotNet;
using static Coil_DHC_Dataset_Maker.dhc_dataset_maker;

namespace Coil_DHC_Dataset_Maker
{
    public static class Program
    {
        // xtodo: 1. cluster strands to see class overlaps.
        // xtodo: 2. measure min,max,average distance between dimorphics/betastrands and LARKs
        // xtodo: 3. 

        public class subseq_match
        {
            public string pdb_id;
            public int chain_number;
            public int alphabetical_chain_number;
            //public List<Atom> pdb_model_chain_atoms;
            public char chain_id;
            public int pdb_model_index;
            public string seq;
            public string subsequence;
            public string dssp;
            public string stride;
            public int total_matches;
            public int array_index;
            public int residue_index;

            //public List<(int array_index, int residue_index)> subseq_indexes;
        }


        public static List<subseq_match> find_subsequences(string pdb_id, string subsequence, bool fix)
        {
            var pdb = Atom.load_atoms_pdb(pdb_id, true, false, false, false, true, true);

            var chain_order = pdb.First().pdb_model_chain_atoms.First().chain_order_in_pdb_file;

            var alphabetical_chain_order = chain_order.OrderBy(a => a.chain_id).Select((a, i) => (array_index: i, chain_id: a.chain_id)).ToList();

            var result = new List<subseq_match>();

            for (var k = 0; k < pdb.Count; k++)
            {
                var master_atoms = Atom.select_amino_acid_master_atoms_ignore_chain(pdb_id, pdb[k].pdb_model_chain_atoms);

                pdb[k].pdb_model_chain_atoms.Clear();
                pdb[k].pdb_model_chain_atoms.AddRange(master_atoms);

            }

            if (fix)
            {

                var seqs = pdb.Select(ch =>
                {
                    var atoms = ch.pdb_model_chain_atoms;
                    var master_atoms = Atom.select_amino_acid_master_atoms_ignore_chain(pdb_id, ch.pdb_model_chain_atoms);
                    var seq = string.Concat(master_atoms.Select(a => a.amino_acid).ToList());
                    return seq;
                }).ToList();


                for (var i = 0; i < seqs.Count; i++)
                {
                    for (var j = 0; j < seqs.Count; j++)
                    {
                        if (i == j) continue;

                        var seqi = seqs[i];
                        var seqj = seqs[j];

                        var chain_id = pdb[i].chain_id;

                        var residsj = pdb[j].pdb_model_chain_atoms.Select(a => a.residue_index).Distinct().OrderBy(a => a).ToList();
                        var residsi = pdb[i].pdb_model_chain_atoms.Select(a => a.residue_index).Distinct().OrderBy(a => a).ToList();

                        // algorithm 1 to deal with missing terminals on one chain
                        if (seqj.Length > seqi.Length && seqj.Contains(seqi))
                        {
                            var pre_final_index = seqj.IndexOf(seqi);
                            var post_first_index = pre_final_index + seqi.Length;

                            var pre = seqj.Substring(0, pre_final_index);
                            var post = seqj.Substring(post_first_index);

                            var pre_len = pre.Length;
                            var post_len = post.Length;



                            var resids_before = residsj.Take(pre_len).OrderBy(a => a).ToList();
                            residsj.Reverse();
                            var resids_after = residsj.Take(post_len).OrderBy(a => a).ToList();
                            residsj.Reverse();

                            var i_first_res_id = residsi.Min();
                            var i_last_res_id = residsi.Max();

                            if (resids_before.Count > 0)
                            {
                                var resids_before_last_id = resids_before.Max();

                            }

                            if (resids_after.Count > 0)
                            {
                                var resids_after_last_id = resids_after.Min();
                            }

                            pdb[i].pdb_model_chain_atoms.InsertRange(0, pdb[j].pdb_model_chain_atoms.Where(a => resids_before.Contains(a.residue_index)).ToList());
                            pdb[i].pdb_model_chain_atoms.AddRange(pdb[j].pdb_model_chain_atoms.Where(a => resids_after.Contains(a.residue_index)).ToList());

                            residsi = pdb[i].pdb_model_chain_atoms.Select(a => a.residue_index).Distinct().OrderBy(a => a).ToList();
                            // note: chain id will still be wrong


                            //Console.WriteLine();
                        }

                        // algorithm 2 to deal with missing termins on one chain
                        var x = residsi.Intersect(residsj).ToList();

                        var pct_intersect = (double)x.Count / (double)((residsi.Count + residsj.Count) / 2);

                        if (pct_intersect >= 0.75)
                        {
                            var i_missing_resids = residsj.Except(residsi).ToList();

                            pdb[i].pdb_model_chain_atoms.AddRange(pdb[j].pdb_model_chain_atoms.Where(a => i_missing_resids.Contains(a.residue_index)).ToList());

                            var z = pdb[i].pdb_model_chain_atoms.OrderBy(a => a.residue_index).ToList();

                            pdb[i].pdb_model_chain_atoms.Clear();
                            pdb[i].pdb_model_chain_atoms.AddRange(z);

                        }
                    }
                }
            }

            for (var k = 0; k < pdb.Count; k++)
            {
                var ch = pdb[k];

                var atoms = ch.pdb_model_chain_atoms;
                var master_atoms = Atom.select_amino_acid_master_atoms_ignore_chain(pdb_id, ch.pdb_model_chain_atoms);

                var seq = string.Concat(master_atoms.Select(a => a.amino_acid).ToList());
                var dssp = string.Concat(master_atoms.Select(a => a.multimer_dssp3).ToList());
                var stride = string.Concat(master_atoms.Select(a => a.multimer_stride3).ToList());

                var match_indexes = seq.IndexOfAll(subsequence);

                if (match_indexes == null || match_indexes.Count == 0) continue;

                foreach (var match_index in match_indexes)
                {
                    var match_atoms = master_atoms.Where((a, i) => i >= match_index && i < match_index + subsequence.Length).ToList();

                    //var match_seq = seq.Substring(match_index, subsequence.Length);
                    //var match_dssp = dssp.Substring(match_index, subsequence.Length);
                    //var match_stride = stride.Substring(match_index, subsequence.Length);

                    var match_seq = string.Concat(match_atoms.Select(a => a.amino_acid).ToList());
                    var match_dssp = string.Concat(match_atoms.Select(a => a.multimer_dssp3).ToList());
                    var match_stride = string.Concat(match_atoms.Select(a => a.multimer_stride3).ToList());

                    var res_id_first = match_atoms.Min(a => a.residue_index);


                    if (match_seq != subsequence) throw new Exception();

                    var match = new subseq_match()
                    {
                        pdb_id = ch.pdb_id,
                        chain_id = ch.chain_id,
                        pdb_model_index = ch.pdb_model_index,
                        chain_number = chain_order.First(a => a.chain_id == ch.chain_id).array_index,
                        alphabetical_chain_number = alphabetical_chain_order.First(a => a.chain_id == ch.chain_id).array_index,
                        seq = seq,
                        subsequence = match_seq,
                        dssp = match_dssp,
                        stride = match_stride,
                        residue_index = res_id_first,
                        array_index = match_index,
                        total_matches = match_indexes.Count

                    };

                    result.Add(match);
                }

            }


            return result;

        }


        public class class_file_line
        {
            public string pair_id;
            public string pdb_id;
            public string dimer_type;
            public string class_name;
            public string symmetry;
            public string parallelism;

            public string strand_seq_merged;
            public string strand_seq_unmerged;
            public string res_id;
            public string chain_colour;

            public string chain_id;

        }

        public static void fill_missing_chains()
        {

            var lines_all = File.ReadAllLines(@"c:\bioinf\fill_missing_chain.csv");

            var lines_header = lines_all.First();

            var lines = lines_all.Skip(1).Select(b =>
            {

                var a = b.Split(',');

                var i = 0;
                return new class_file_line()
                {
                    pair_id = a[i++].Trim(),
                    pdb_id = a[i++].Trim(),
                    dimer_type = a[i++].Trim(),
                    class_name = a[i++].Trim(),
                    symmetry = a[i++].Trim(),
                    parallelism = a[i++].Trim(),
                    strand_seq_merged = a[i++].Trim(),
                    strand_seq_unmerged = a[i++].Trim(),
                    res_id = a[i++].Trim(),
                    chain_colour = a[i++].Trim()

                };

            }).Distinct().ToList();

            //0 PDB ID
            //1 Dimer Type
            //2 Class
            //3 Symmetry
            //4 Parallelism
            //5 Strand Sequence
            //6 Res ID

            // first, fill in the 'missing' values, since those are the 'first' 
            // that won't work, because it might be a heterodimer, with different sequences, but both containing the strand subsequence

            // 

            //lines = lines.Skip(lines.FindIndex(a=>a.pdb_id.Equals("2NS9", StringComparison.CurrentCultureIgnoreCase))).ToList();
            //lines = lines.Where(a => a.res_id.Length > 0).ToList();


            var errors = 0;

            var tasks = new List<Task>();

            for (var i = 0; i < lines.Count; i++)
            {
                var line = lines[i];


                Console.WriteLine(line.pdb_id);

                var task = Task.Run(() =>
                {
                    var fix = false;

                    do
                    {


                        var lookup = find_subsequences(line.pdb_id, line.strand_seq_merged, fix);
                        fix = false;

                        // filter lookup results by res id, if specified
                        if (line.res_id.Length > 0)
                        {

                            var lookup2 = lookup.Where(a => a.residue_index == int.Parse(line.res_id)).ToList();

                            if (lookup2.Count == 0) throw new Exception();


                            lookup = lookup2;
                        }


                        // filter by pymol chain colour, 0 is red, 1 is pink

                        if (line.chain_colour == "Red" && lookup.Any(a => a.chain_number == 0))
                        {
                            var colour_lookup = lookup.Where(a => a.chain_number == 0).ToList();

                            if (colour_lookup.Count > 1)
                                Console.WriteLine();

                            line.res_id = colour_lookup.First().residue_index.ToString();
                            line.chain_id = colour_lookup.First().chain_id.ToString();
                        }

                        else if (line.chain_colour == "Pink" && lookup.Any(a => a.chain_number == 1))
                        {
                            var colour_lookup = lookup.Where(a => a.chain_number == 1).ToList();

                            if (colour_lookup.Count > 1)
                                Console.WriteLine();

                            line.res_id = colour_lookup.First().residue_index.ToString();
                            line.chain_id = colour_lookup.First().chain_id.ToString();
                        }

                        else
                        {
                            errors++;
                            Console.WriteLine("errors: " + errors + ", fix: " + fix);
                            //throw new Exception();
                            fix = true;
                        }

                    } while (fix);

                });

                tasks.Add(task);
            }

            Task.WaitAll(tasks.ToArray());

            var result = lines.Select(a => string.Join(",", new string[] { a.pair_id, a.pdb_id, a.dimer_type, a.class_name, a.symmetry, a.parallelism, a.strand_seq_merged, a.strand_seq_unmerged, a.res_id, a.chain_colour, a.chain_id })).ToList();
            result.Insert(0, lines_header + "," + "chain_id");

            var fn = @"c:\bioinf\fixed_dataset.csv";

            Directory.CreateDirectory(Path.GetDirectoryName(fn));

            File.WriteAllLines(fn, result);
            Console.WriteLine("finished.");
            Console.ReadLine();

        }

        //public static void get_missing_data()
        //{
        //    var lines = File.ReadAllLines(@"c:\bioinf\dimorphics_data_set.csv").Skip(1).Select(a => a.Split(',').ToList()).ToList();

        //    var pdb_id_col = 1;
        //    var seq_col = 5;

        //    var list = lines.Select(a => (a[pdb_id_col], a[seq_col], a)).ToList();

        //    //var results = new List<string>();

        //    var tasks = new List<Task<List<string>>>();

        //    var max_tasks = 1000;

        //    foreach (var seq in list)
        //    {
        //        var z = seq;

        //        var task = Task.Run(() =>
        //        {
        //            var x = find_subsequences(z.Item1, z.Item2);

        //            for (var i = 0; i < x.Count; i++)
        //            {
        //               // x[i] = string.Join(",", z.a) + "," + x[i];
        //            }
        //            //results.AddRange(x);

        //            return x;

        //        });

        //        tasks.Add(task);

        //        while (tasks.Count(a => !a.IsCompleted) >= max_tasks)
        //        {
        //            Task.WaitAny(tasks.ToArray<Task>());
        //        }
        //    }

        //    Task.WaitAll(tasks.ToArray<Task>());


        //    var results = tasks.SelectMany(a => a.Result).ToList();

        //    File.WriteAllLines(@"c:\bioinf\res_ids_list.csv", results);

        //    Console.WriteLine("Finished.");
        //    Console.ReadLine();
        //}

        // xtodo: encode the foldx scores
        internal static object _console_lock = new object();

        public static void Write(string text = "")
        {
            lock (_console_lock)
            {
                Console.Write(text);
            }
        }

        public static void WriteLine(string text = "")
        {
            lock (_console_lock)
            {
                Console.WriteLine(text);
            }
        }

        //subsequence_classification_data.encode_sequence("x", "ALGY");

        //xtodo: create dhc dataset with dhc found by dssp3, dssp7, stride3, stride7, for comparison
        //xtodo: to see which secondary structure definition gives the best performance in the classifier


        //xtodo: ensure that data is placed in the correct column, in case of any features being empty/null/not included in some rows.

        //var pdb_files = Directory.GetFiles(@"c:\betastrands_dataset\pdb\", "*.pdb");
        //pdb_files.ToList().ForEach(a=>Atom.run_dssp(Path.GetFileNameWithoutExtension(a)));
        //pdb_files.ToList().ForEach(a=>Atom.split_pdb_chains(Path.GetFileNameWithoutExtension(a)));
        //pdb_files.ToList().ForEach(a=>Program.WriteLine($"dssp -i {Path.GetFileNameWithoutExtension(a)}.pdb -o {Path.GetFileNameWithoutExtension(a)}.dssp"));
        //return;

        //if (run_larks_unknown)
        //{
        //    class_list.Add(
        //        (
        //        nameof(run_larks_unknown).EndsWith("unknown") ? +1 : -1,
        //        nameof(run_larks_unknown).Replace("run_", ""),
        //        $@"c:\bioinf\e_{nameof(run_larks_unknown).Replace("run_", "")}{(sequence_only_features ? "_seq_only" : "")}_{uid}.txt",
        //        $@"c:\bioinf\d_{nameof(run_larks_unknown).Replace("run_", "")}{(sequence_only_features ? "_seq_only" : "")}_{uid}.txt")

        //    );
        //}

        //else if (class_name == nameof(run_larks_unknown).Replace("run_", ""))
        //{
        //    var larks = new List<string>() { "STGGYG", "GYNGFG", "SYSGYS", "SYSSYGQS", "GFGNFGTS" };
        //    sequence_list = larks;
        //}

        public class class_info
        {
            public int class_id;
            public substructure_type substructure_type;
            public int min_subsequence_length;
            public string class_name;
            //public string class_output_file;
            //public string class_dm_output_file;
            //public string class_comments_output_file;
        }

        public static void extract_all_pdbs()
        {
            var pdbs = Directory.GetFiles(@"c:\betastrands_dataset\pdb\", "*.pdb").ToList();

            foreach (var p in pdbs)
            {
                Console.WriteLine("Extracting " + p);
                Atom.extract_split_pdb_chains(p, null);
            }
        }

        public static void repair_all_extracted_pdbs(bool run)
        {
            var split_pdbs = Directory.GetFiles(@"c:\betastrands_dataset\pdb_split\", "*.pdb").ToList();

            var tasks = new List<Task>();

            foreach (var pp in split_pdbs)
            {
                var p = pp;

                var task = Task.Run(() =>
                {
                    Console.WriteLine("Repairing " + p);
                    foldx_caller.foldx_repair_pdb(Path.GetFileNameWithoutExtension(p), run, Path.GetDirectoryName(p) + @"\");
                });

                tasks.Add(task);

                while (tasks.Count(a => !a.IsCompleted) >= Environment.ProcessorCount)
                {
                    Task.WaitAny(tasks.ToArray());
                }
            }

            Task.WaitAll(tasks.ToArray());

        }

        public static List<string> get_pdb_sequences(bool limited_to_dimorphics_and_standard = true)
        {

            var dimorphics_data = File.ReadAllLines(@"C:\betastrands_dataset\csv\distinct dimorphics list.csv")
                .Skip(1).Where(a => !string.IsNullOrWhiteSpace(a.Replace(",", ""))).Select((a, i) =>
                {
                    var x = a.Split(',');
                    return (
                        pdb_id: x[0].ToUpperInvariant(),
                        dimer_type: x[1],
                        class_name: x[2],
                        symmetry_mode: x[3],
                        parallelism: x[4],
                        chain_number: int.Parse(x[5]) - 1,
                        strand_seq: x[6],
                        optional_res_index: x[7]
                    );
                }).ToList();

            if (limited_to_dimorphics_and_standard)
            {
                dimorphics_data = dimorphics_data.Where(a => (a.class_name == "Dimorphic" || a.class_name == "Single") || (a.class_name == "Standard" || a.class_name == "Multiple")).ToList();
            }


            var pdb_id_list = dimorphics_data.Select(a => a.pdb_id.Trim().ToUpperInvariant()).Distinct().ToList();


            var pdb_folder = @"c:\betastrands_dataset\pdb\";


            var pdbs = Directory.GetFiles(pdb_folder, "*.pdb").ToList();

            var seqs = new List<string>();

            var reverse_lookup = new List<(string pdb, string seq)>();
            foreach (var pdb_id in pdbs)
            {
                var pdb_id_code = Path.GetFileNameWithoutExtension(pdb_id).Trim().ToUpperInvariant().Substring(0, 4);

                //if (!pdb_id_list.Contains(pdb_id_code)) continue;

                Console.WriteLine("Extracting sequence from " + pdb_id);


                var chains = Atom.load_atoms_pdb(pdb_id, true, false, false, false, false, false, pdb_folder);

                foreach (var chain in chains)
                {


                    var pdb_seqs = Atom.amino_acid_sequence(pdb_id, chain.pdb_model_chain_atoms);

                    foreach (var s in pdb_seqs)
                    {
                        seqs.Add(s.aa_sequence);
                        reverse_lookup.Add((chain.pdb_id + chain.chain_id, s.aa_sequence));
                    }
                }
            }

            var output_filename = @"c:\bioinf\betastrands_dataset_sequences.txt";
            var rename_filename = @"c:\bioinf\rename_num_to_pdb.bat";


            if (limited_to_dimorphics_and_standard)
            {
                output_filename = @"c:\bioinf\betastrands_dataset_sequences_limited.txt";
                rename_filename = @"c:\bioinf\rename_num_to_pdb_limited.bat";

            }


            Directory.CreateDirectory(Path.GetDirectoryName(output_filename));

            seqs = seqs.OrderBy(a => a).Distinct().ToList();

            //            var rename_list = seqs.SelectMany((a, i) => reverse_lookup.Where(b=>b.seq==a).Select(b=> $"copy {i}.* {b.pdb}.*")).ToList();

            var rename_list = reverse_lookup.Select(a => (a.pdb, a.seq, index: seqs.IndexOf(a.seq), cmd: new[] { $"del {a.pdb}.*", $"copy {seqs.IndexOf(a.seq)}.* {a.pdb}.*" })).ToList();

            File.WriteAllLines(output_filename, seqs);
            File.WriteAllLines(rename_filename, rename_list.SelectMany(a => a.cmd).ToList());

            //for (var i = 0; i < seqs.Count; i++)
            //{
            //    var file_seq = new List<string>();

            //    file_seq.Add(">" + i);
            //    file_seq.Add(seqs[i]);

            //    File.WriteAllLines(@"C:\betastrands_dataset\iupred2a\" + i + ".seq", file_seq);
            //}

            return seqs;

        }



        public static void get_pssms(string blast_db = "swissprot", string db_folder = @"c:\blast\db\", bool remote = false, bool run = false)
        {
            int num_iterations = 3;

            // this method:
            // 0. loads 'txt_sequences' text file of dataset sequences (but only uses those sequences also found within 'txt_sequences_limited')
            // 1. if 'output_fasta_files' is true, outputs the fasta sequences in individual files
            // 2. makes a list of command lines to run blast, skipping any pssm files which already exist (& can also run blast if 'run' set to true)

            var txt_sequences = File.ReadAllLines(@"c:\betastrands_dataset\betastrands_dataset_sequences.txt");
            var txt_sequences_limited = File.ReadAllLines(@"c:\betastrands_dataset\betastrands_dataset_sequences_limited.txt");

            //var excluded = txt_sequences.Except(txt_sequences_limited).ToList();
            //var included = txt_sequences.Except(excluded).ToList();

            //var fasta_file = $@"c:\betastrands_dataset\fasta\{index1}.fasta";
            //var pssm_file = $@"c:\betastrands_dataset\blast_pssm_{blast_db}_{num_iterations}_{(remote ? "remote" : "local")}\{index1}.pssm";


            var fasta_blast_input_folder = $@"c:\fasta\";
            var fasta_output_folder = $@"c:\fasta\";

            //var pssm_blast_output_folder = $@"h:\data\blast_pssm_{blast_db}_{(!remote ? num_iterations.ToString() + "_" : "")}{(remote ? "new_remote" : "local")}";
            //var already_calculated_pssms = File.ReadAllLines(@"C:\Users\aaron\Desktop\{}");





            var output_fasta_files = false;
            if (output_fasta_files)
            {
                for (var index = 0; index < txt_sequences.Length; index++)
                {
                    var fasta_file = Path.Combine(fasta_output_folder, $@"{index}.fasta");

                    if (File.Exists(fasta_file) && new FileInfo(fasta_file).Length > 0) continue;

                    Console.WriteLine($"Saving fasta {index}");
                    var t_seq = txt_sequences[index];
                    var seq = new Sequence($">seq_{index}", t_seq);
                    Sequence.Save(fasta_file, seq);
                }
            }

            //var tasks = new List<Task<List<string>>>();
            //var max_tasks = 1; // Environment.ProcessorCount
            var skipped_not_limited = 0;
            var skipped_pssm_exists = 0;

            var dna_pred = new List<string>();

            var emax = 100.0;
            var emin = 0.00001;
            var ediv = 10;

            for (var evalue = emax; evalue >= emin; evalue /= ediv)
            {
                for (var inclusion_ethresh = emax; inclusion_ethresh >= emin; inclusion_ethresh /= ediv)
                {
                    var cmd_lines = new List<string>() { };

                    var pssm_blast_input_folder = $@"c:\pssm\blast_pssm_{(!remote ? num_iterations.ToString() + "_" : "")}{(remote ? "remote" : "local")}_{blast_db}_{evalue.ToString("000.00000", CultureInfo.InvariantCulture)}_{inclusion_ethresh.ToString("000.00000", CultureInfo.InvariantCulture)}";
                    var pssm_blast_output_folder = $@"c:\pssm\blast_pssm_{(!remote ? num_iterations.ToString() + "_" : "")}{(remote ? "remote" : "local")}_{blast_db}_{evalue.ToString("000.00000", CultureInfo.InvariantCulture)}_{inclusion_ethresh.ToString("000.00000", CultureInfo.InvariantCulture)}";

                    Directory.CreateDirectory(pssm_blast_input_folder);
                    Directory.CreateDirectory(pssm_blast_output_folder);

                    for (var index = 0; index < txt_sequences.Length; index++)
                    {
                        var index1 = index;

                        if (!txt_sequences_limited.Contains(txt_sequences[index1]))
                        {
                            skipped_not_limited++;
                            continue;
                        }
                        //Console.WriteLine($"Querying psi blast with sequence {index1}");

                        var dna_str = $@"python3 predict.py {index1} {txt_sequences[index1]}";
                        if (!dna_pred.Contains(dna_str)) dna_pred.Add(dna_str);

                        var fasta_file = Path.Combine(fasta_blast_input_folder, $@"{index1}.fasta");

                        var pssm_local_file = Path.Combine(pssm_blast_input_folder, $@"{index1}.pssm");

                        var pssm_copy_file = Path.Combine(pssm_blast_output_folder, $@"{index1}.pssm");

                        if (File.Exists(pssm_local_file) && new FileInfo(pssm_local_file).Length > 0 && File.ReadAllLines(pssm_local_file).Where(a => !string.IsNullOrWhiteSpace(a)).ToList().Any(x => x.StartsWith("PSI Gapped")))
                        {
                            skipped_pssm_exists++;
                            continue;
                        }



                        var cmd_line = pssm.run_psi_blast_get_pssm(fasta_file, pssm_local_file, blast_db, db_folder, num_iterations, evalue.ToString("000.00000", CultureInfo.InvariantCulture), inclusion_ethresh.ToString("000.00000", CultureInfo.InvariantCulture), remote, run);

                        cmd_lines.Add(cmd_line);

                        if (!string.Equals(pssm_local_file, pssm_copy_file, StringComparison.Ordinal))
                        {
                            cmd_lines.Add($"copy /y \"{pssm_local_file}\" \"{pssm_copy_file}\"");
                        }



                    }

                    var batch_pssm_local_file = Path.Combine(pssm_blast_input_folder, $@"blast_pssm_{(!remote ? num_iterations.ToString() + "_" : "")}{(remote ? "remote" : "local")}_{blast_db}_{evalue.ToString("000.00000", CultureInfo.InvariantCulture)}_{inclusion_ethresh.ToString("000.00000", CultureInfo.InvariantCulture)}.bat");
                    File.WriteAllLines(batch_pssm_local_file, cmd_lines);

                }
            }

            //Task.WaitAll(tasks.ToArray());

            //var all_cmd_lines = tasks.SelectMany(a => a.Result).ToList();

            File.WriteAllLines(@"c:\betastrands_dataset\dna_pred.bat", dna_pred);

            Console.WriteLine($@"{blast_db}_{num_iterations}_{(remote ? "remote" : "local")}: {nameof(skipped_not_limited)}:{skipped_not_limited} {nameof(skipped_pssm_exists)}:{skipped_pssm_exists}");
        }

        public static void load_pssms()
        {
            var txt_sequences = File.ReadAllLines(@"c:\bioinf\betastrands_dataset_sequences.txt");



            for (var index = 0; index < txt_sequences.Length; index++)
            {
                var index1 = index;

                Console.WriteLine($"Loading pssm {index1}");


                var pssm_matrix = pssm.load_psi_blast_pssm($@"c:\bs_pssm\{index1}.pssm");

                pssm_matrix = pssm.normalise_pssm(pssm_matrix);

                Console.WriteLine($"Loaded pssm {index1}");

            }


        }

        public static (string header, string sequence) get_uniprot_sequence(string uniprot_id)
        {
            var uniprot_url = $@"https://www.uniprot.org/uniprot/{uniprot_id}.fasta";


            var wc = new WebClient();
            var fasta = wc.DownloadString(uniprot_url).Split(new char[] { '\r', '\n' }, StringSplitOptions.RemoveEmptyEntries);

            var header = string.Join("", fasta.Where(b => b.StartsWith(">")).ToList());
            var sequence = string.Join("", fasta.Where(b => !b.StartsWith(">")).ToList());

            //Console.WriteLine(sequence);
            return (header, sequence);
        }

        public static void get_uniprot_sequences(bool limited_to_dimorphics_and_standard = true)
        {
            //var dimorphics_data = File.ReadAllLines(@"C:\betastrands_dataset\csv\distinct dimorphics list.csv")
            //    .Skip(1).Where(a => !string.IsNullOrWhiteSpace(a.Replace(",", ""))).Select((a, i) =>
            //    {
            //        var x = a.Split(',');
            //        return (
            //            pdb_id: x[0].ToUpperInvariant(),
            //            dimer_type: x[1],
            //            class_name: x[2],
            //            symmetry_mode: x[3],
            //            parallelism: x[4],
            //            chain_number: int.Parse(x[5]) - 1,
            //            strand_seq: x[6],
            //            optional_res_index: x[7]
            //        );
            //    }).ToList();

            //if (limited_to_dimorphics_and_standard)
            //{
            //    dimorphics_data = dimorphics_data.Where(a => a.class_name == "Single" || a.class_name == "Multiple").ToList();
            //}


            //var pdb_id_list = dimorphics_data.Select(a => a.pdb_id.Trim().ToUpperInvariant()).Distinct().ToList();


            var sifts_lookup = File.ReadAllLines(@"c:\betastrands_dataset\sifts\pdb_chain_uniprot.csv").Skip(2).Select(a =>
            {

                var x = a.Split(',');
                var i = 0;
                var y = (
                    PDB: x[i++].Substring(0, 4).ToUpperInvariant(),
                    CHAIN: x[i++][0],
                    SP_PRIMARY: x[i++],
                    RES_BEG: x[i++],
                    RES_END: x[i++],
                    PDB_BEG: x[i++],
                    PDB_END: x[i++],
                    SP_BEG: x[i++],
                    SP_END: x[i++]);

                return y;
            }).ToList();

            var pdb_folder = @"c:\betastrands_dataset\pdb_split\";
            var uniprot_folder = @"c:\betastrands_dataset\uniprot\";


            var files = Directory.GetFiles(pdb_folder, "*.pdb");

            foreach (var a in files)
            {
                Console.WriteLine(a);

                var file = Path.GetFileNameWithoutExtension(a).ToUpperInvariant();
                var pdb = file.Substring(0, 4);
                var chain = file[4];
                var uniprot_id = sifts_lookup.FirstOrDefault(b => b.PDB.Equals(pdb, StringComparison.InvariantCultureIgnoreCase) && b.CHAIN == chain).SP_PRIMARY;

                if (string.IsNullOrWhiteSpace(uniprot_id)) continue;

                var uniprot_fasta = get_uniprot_sequence(uniprot_id);

                //var data = (pdb: pdb, chain: chain, uniprot_id: uniprot_id, uniprot_header:uniprot_fasta.header, uniprot_sequence:uniprot_fasta.sequence);

                Directory.CreateDirectory(uniprot_folder);
                File.WriteAllLines(Path.Combine(uniprot_folder, pdb + chain + ".fasta"), new string[] { uniprot_fasta.header, uniprot_fasta.sequence });
            }






        }

        public static void check_buildmodel_position_scan_files_are_complete()
        {
            var cmd_list_file = @"C:\betastrands_dataset\foldx\bat\foldx_calc_buildmodel_position_scan.bat";
            var cmd_list_file_missing = @"C:\betastrands_dataset\foldx\bat\missing_foldx_calc_buildmodel_position_scan.bat";
            var cmd_list_file_missing_del = @"C:\betastrands_dataset\foldx\bat\del_missing_foldx_calc_buildmodel_position_scan.bat";

            var cmd_list = File.ReadAllLines(cmd_list_file).Where(a => a.Contains("--mutant-file=")).ToList();


            var cmd_list_parsed = cmd_list.Select(a =>
            {
                var s = a.Split();


                var mutant_file = s.First(b => b.StartsWith("--mutant-file=")).Split('=').Last();
                if (!File.Exists(mutant_file) || new FileInfo(mutant_file).Length == 0)
                {
                    Console.WriteLine($"missing mutant file: {mutant_file}");
                }

                var dif_fxout_file = s[3].Replace("\"", "");
                var average_fxout_file = s[3].Replace("\"", "").Replace("Dif", "Average");
                var raw_fxout_file = s[3].Replace("\"", "").Replace("Dif", "Raw");



                var mutant_data = File.ReadAllLines(mutant_file).Where(b => !string.IsNullOrWhiteSpace(b)).Select(b => b.Trim()).ToList();


                var dif_fxout_data = (!File.Exists(dif_fxout_file) || new FileInfo(dif_fxout_file).Length == 0) ? new List<string>() : File.ReadAllLines(dif_fxout_file).Skip(9).Where(b => !string.IsNullOrWhiteSpace(b)).ToList();
                var average_fxout_data = (!File.Exists(average_fxout_file) || new FileInfo(average_fxout_file).Length == 0) ? new List<string>() : File.ReadAllLines(average_fxout_file).Skip(9).Where(b => !string.IsNullOrWhiteSpace(b)).ToList();
                var raw_fxout_data = (!File.Exists(raw_fxout_file) || new FileInfo(raw_fxout_file).Length == 0) ? new List<string>() : File.ReadAllLines(raw_fxout_file).Skip(9).Where(b => !string.IsNullOrWhiteSpace(b)).ToList();


                return (dif_fxout_file, dif_fxout_data, average_fxout_file, average_fxout_data, raw_fxout_file, raw_fxout_data, mutant_file, mutant_data, a);

            }).ToList();

            var missing = cmd_list_parsed.Where(a => a.mutant_data.Count != a.dif_fxout_data.Count).ToList();
            var not_missing = cmd_list_parsed.Except(missing).ToList();

            foreach (var m in missing)
            {
                var total_missing = m.mutant_data.Count - m.dif_fxout_data.Count;
                Console.WriteLine($"missing {total_missing} / {m.mutant_data.Count}: {m.mutant_file}");
            }



            Console.WriteLine($"total missing: {missing.Count}");
            Console.WriteLine($"total not missing: {not_missing.Count}");






            var missing_list = missing.Select(a => a.a).ToList();
            var missing_list_del = missing_list.Where(a => a.StartsWith("if not exist", StringComparison.InvariantCultureIgnoreCase)).Select(b => "del " + b.Split()[3]).ToList();
            missing_list_del.AddRange(missing_list_del.Select(a => a.Replace("Dif", "Average")).ToList());
            missing_list_del.AddRange(missing_list_del.Select(a => a.Replace("Dif", "Raw")).ToList());

            File.WriteAllLines(cmd_list_file_missing, missing_list);
            File.WriteAllLines(cmd_list_file_missing_del, missing_list_del);


            Console.WriteLine();

        }


        public static void check_buildmodel_subseq_mutant_files_are_complete()
        {
            var cmd_list_file = @"C:\betastrands_dataset\foldx\bat\foldx_calc_buildmodel_subsequence_mutant.bat";
            var cmd_list_file_missing = @"C:\betastrands_dataset\foldx\bat\missing_foldx_calc_buildmodel_subsequence_mutant.bat";
            var cmd_list_file_missing_del = @"C:\betastrands_dataset\foldx\bat\del_missing_foldx_calc_buildmodel_subsequence_mutant.bat";

            var cmd_list = File.ReadAllLines(cmd_list_file).Where(a => a.Contains("--mutant-file=")).ToList();



            var cmd_list_parsed = cmd_list.Select(a =>
            {
                var s = a.Split();


                var mutant_file = s.First(b => b.StartsWith("--mutant-file=")).Split('=').Last();
                if (!File.Exists(mutant_file) || new FileInfo(mutant_file).Length == 0)
                {
                    Console.WriteLine($"missing mutant file: {mutant_file}");
                }

                var dif_fxout_file = s[3].Replace("\"", "");
                var average_fxout_file = s[3].Replace("\"", "").Replace("Dif", "Average");
                var raw_fxout_file = s[3].Replace("\"", "").Replace("Dif", "Raw");



                var mutant_data = File.ReadAllLines(mutant_file).Where(b => !string.IsNullOrWhiteSpace(b)).Select(b => b.Trim()).ToList();


                var dif_fxout_data = (!File.Exists(dif_fxout_file) || new FileInfo(dif_fxout_file).Length == 0) ? new List<string>() : File.ReadAllLines(dif_fxout_file).Skip(9).Where(b => !string.IsNullOrWhiteSpace(b)).ToList();
                var average_fxout_data = (!File.Exists(average_fxout_file) || new FileInfo(average_fxout_file).Length == 0) ? new List<string>() : File.ReadAllLines(average_fxout_file).Skip(9).Where(b => !string.IsNullOrWhiteSpace(b)).ToList();
                var raw_fxout_data = (!File.Exists(raw_fxout_file) || new FileInfo(raw_fxout_file).Length == 0) ? new List<string>() : File.ReadAllLines(raw_fxout_file).Skip(9).Where(b => !string.IsNullOrWhiteSpace(b)).ToList();


                return (dif_fxout_file, dif_fxout_data, average_fxout_file, average_fxout_data, raw_fxout_file, raw_fxout_data, mutant_file, mutant_data, a);

            }).ToList();

            var missing = cmd_list_parsed.Where(a => a.mutant_data.Count != a.dif_fxout_data.Count).ToList();
            var not_missing = cmd_list_parsed.Except(missing).ToList();

            foreach (var m in missing)
            {
                var total_missing = m.mutant_data.Count - m.dif_fxout_data.Count;
                Console.WriteLine($"missing {total_missing} / {m.mutant_data.Count}: {m.mutant_file}");
            }



            Console.WriteLine($"total missing: {missing.Count}");
            Console.WriteLine($"total not missing: {not_missing.Count}");


            var missing_list = missing.Select(a => a.a).ToList();
            var missing_list_del = missing_list.Where(a => a.StartsWith("if not exist", StringComparison.InvariantCultureIgnoreCase)).Select(b => "del " + b.Split()[3]).ToList();
            missing_list_del.AddRange(missing_list_del.Select(a => a.Replace("Dif", "Average")).ToList());
            missing_list_del.AddRange(missing_list_del.Select(a => a.Replace("Dif", "Raw")).ToList());

            File.WriteAllLines(cmd_list_file_missing, missing_list);
            File.WriteAllLines(cmd_list_file_missing_del, missing_list_del);


            Console.WriteLine();

        }

        public static void check_ala_scanning_files_are_complete()
        {
            var ala_scan_folder = @"C:\betastrands_dataset\foldx_ala_scan\";

            var cmd_list_file = @"C:\betastrands_dataset\foldx\bat\foldx_calc_ala_scanning.bat.skip";
            var cmd_list_file_missing = @"C:\betastrands_dataset\foldx\bat\missing_foldx_calc_ala_scanning.bat.skip";
            var cmd_list_file_missing_del = @"C:\betastrands_dataset\foldx\bat\del_missing_foldx_calc_ala_scanning.bat.skip";

            var cmd_list = File.ReadAllLines(cmd_list_file).Where(a => a.Contains("--pdb=")).ToList();

            var cmd_list_parsed = cmd_list.Select(a =>
            {
                var s = a.Split();

                var id = Path.GetFileNameWithoutExtension(s[3].Replace("\"", "").Replace("--pdb=", ""));

                var pdb_file = @"C:\betastrands_dataset\pdb_split_repair\" + id + ".pdb";

                var ala_file = Path.Combine(ala_scan_folder, id + "_AS.fxout");

                var pdb_res = File.ReadAllLines(pdb_file).Where(b => b.StartsWith("ATOM")).Select(c =>
                    c.Substring(17, 3) + " " + int.Parse(c.Substring(22, 4).Trim())).Distinct().ToList();

                var ala_res = File.ReadAllLines(ala_file).Where(c => !string.IsNullOrWhiteSpace(c))
                    .Select(c => c.Substring(0, c.IndexOf(" to ")).Trim()).ToList();

                ala_res = ala_res.Select(b => b.Replace("H1S", "HIS")).ToList();
                ala_res = ala_res.Select(b => b.Replace("H2S", "HIS")).ToList();

                var present = pdb_res.Intersect(ala_res).ToList();
                var missing1 = pdb_res.Except(ala_res).ToList();
                var extra = ala_res.Except(pdb_res).ToList();

                return (id: id, pdb_res: pdb_res, ala_res: ala_res, present: present, missing1: missing1, extra: extra, a);



                //positions: s[5].Split('=').Last().Split(','));
            }).ToList();

            var missing = cmd_list_parsed.Where(a => a.missing1.Count > 0 || a.extra.Count > 0).ToList();

            Console.WriteLine("missed: " + missing.Count);
            foreach (var m in missing)
            {
                Console.WriteLine($@"{m.id}, extra: {string.Join("|", m.extra)}, missing: {string.Join("|", m.missing1)}");
            }


            var missing_list = missing.Select(a => a.a).ToList();
            var missing_list_del = missing_list.Where(a => a.StartsWith("if not exist", StringComparison.InvariantCultureIgnoreCase)).Select(b => "del " + b.Split()[3]).ToList();
            File.WriteAllLines(cmd_list_file_missing, missing_list);
            File.WriteAllLines(cmd_list_file_missing_del, missing_list_del);

            Console.WriteLine();
        }

        public static void check_foldx_position_scan_files_are_completed()
        {

            var missing_res = compare_repaired_pdb_to_pdb();

            var ps_folder = @"C:\betastrands_dataset\foldx\ps\";
            var file_mask = "PS_*.txt";

            var ps_files = Directory.GetFiles(ps_folder, file_mask, SearchOption.TopDirectoryOnly);

            var cmd_list_file = @"C:\betastrands_dataset\foldx\bat\foldx_calc_position_scanning.bat";
            var cmd_list_file_missing = @"C:\betastrands_dataset\foldx\bat\missing_foldx_calc_position_scanning.bat";
            var cmd_list_file_missing_del = @"C:\betastrands_dataset\foldx\bat\del_missing_foldx_calc_position_scanning.bat";

            var cmd_list = File.ReadAllLines(cmd_list_file).Where(a => a.Contains("--positions=")).ToList();

            var cmd_list_parsed = cmd_list.Select(a =>
            {
                var s = a.Split();

                return (file: s[3].Replace("\"", ""), positions: s[5].Split('=').Last().Split(','), a);
            }).ToList();

            var incomplete_files = new List<string>();
            var complete_files = new List<string>();

            var missing = new List<string>();

            foreach (var c in cmd_list_parsed)
            {
                var file = c.file;

                var pdb = Path.GetFileNameWithoutExtension(file).Substring(3, 5).Trim();

                var pdb_missing_res = missing_res.Where(a => string.Equals(a.id, pdb, StringComparison.InvariantCultureIgnoreCase)).ToList();

                var actual_miss = pdb_missing_res.SelectMany(a => a.original_except_repair).ToList();

                var x = c.positions.Select(a => a[0] + a.Substring(1)).ToList();

                // convert 'd' to the 25 codes

                x = x.SelectMany(a =>
                {
                    var rid = int.Parse(string.Join("", a.Where(b => char.IsDigit(b)).ToList()));

                    if (actual_miss.Contains(rid))
                    {

                        return new List<string>() { };
                    }

                    if (a.Last() == 'd')
                    {                //123456789012345678901234567890
                        var d_codes = "GALVIPRTSCMKEQDNWYFypsh"; //e

                        var output = d_codes.Select(d => a.Substring(0, a.Length - 1) + d).ToList();
                        return output;
                    }
                    else
                    {
                        return new List<string>() { a };

                    }
                }).Where(a => !string.IsNullOrWhiteSpace(a)).Distinct().ToList();

                var x_res = x.Select(a => int.Parse(string.Join("", a.Where(b => char.IsDigit(b)).ToList()))).Distinct().ToList();

                var scanning_data = File.Exists(c.file) && new FileInfo(c.file).Length > 0 ? File.ReadAllLines(c.file).Select(a =>
                    {
                        var z = a.Split();
                        z[0] = foldx_caller.foldx_residues_aa_mutable.First(d => d.foldx_aa_code3 == z[0].Substring(0, 3)).standard_aa_code1 + z[0].Substring(3);
                        return z;
                    }).Distinct().ToList() : new List<string[]>();


                var data_codes = scanning_data.Select(a => a.First()).Distinct().ToList();
                var data_codes_res = data_codes.Select(a => int.Parse(string.Join("", a.Where(b => char.IsDigit(b)).ToList()))).Distinct().ToList();

                var missing1 = x.Except(data_codes).Distinct().ToList();

                if (missing1.Count > 0)
                {
                    missing.Add(c.a);

                    incomplete_files.Add(file);

                    Console.WriteLine();
                    Console.WriteLine("File: " + file);
                    Console.WriteLine();
                    Console.WriteLine("Missing " + missing1.Count + ": " + string.Join(", ", missing1));
                    Console.WriteLine();
                }
                else
                {
                    complete_files.Add(file);

                }

            }

            Console.WriteLine("incomplete_file_count: " + incomplete_files.Count);
            Console.WriteLine("complete_file_count: " + complete_files.Count);

            //File.WriteAllLines(@"C:\betastrands_dataset\foldx\ps_incomplete_files.txt", incomplete_files);



            var missing_list = missing.Select(a => a).ToList();
            var missing_list_del = missing_list.Where(a => a.StartsWith("if not exist", StringComparison.InvariantCultureIgnoreCase)).Select(b => "del " + b.Split()[3]).ToList();
            File.WriteAllLines(cmd_list_file_missing, missing_list);
            File.WriteAllLines(cmd_list_file_missing_del, missing_list_del);

            Console.WriteLine();
        }

        public static List<(string id, string id_repair, List<int> res_id_original, List<int> res_id_repair, List<int> intersect, List<int> original_except_repair, List<int> repair_except_original)> compare_repaired_pdb_to_pdb()
        {
            var split_pdbs = @"C:\betastrands_dataset\pdb_split\";
            var split_pdbs_repair = @"C:\betastrands_dataset\pdb_split_repair\";

            var pdbs_original = Directory.GetFiles(split_pdbs, "*.pdb");
            var pdbs_repair = Directory.GetFiles(split_pdbs_repair, "*.pdb");

            var ids = pdbs_original.Select(a => Path.GetFileNameWithoutExtension(a)).ToList();

            var pdbs = ids.Select((a, i) => (id: Path.GetFileNameWithoutExtension(pdbs_original[i]), id_repair: Path.GetFileNameWithoutExtension(pdbs_repair[i]), pdb_original: pdbs_original[i], pdb_repair: pdbs_repair[i])).ToList();

            List<(string id, string id_repair, List<int> res_id_original, List<int> res_id_repair, List<int> intersect, List<int> original_except_repair, List<int> repair_except_original)> diff = pdbs.Select(a =>
            {
                var o = File.ReadAllLines(a.pdb_original).Where(b => b.StartsWith("ATOM")).Select(b => int.Parse(b.Substring(22, 4).Trim())).Distinct().ToList();
                var r = File.ReadAllLines(a.pdb_repair).Where(b => b.StartsWith("ATOM")).Select(b => int.Parse(b.Substring(22, 4).Trim())).Distinct().ToList();
                var i = o.Intersect(r).ToList();
                var e1 = o.Except(r).ToList();
                var e2 = r.Except(o).ToList();

                return (id: a.id, id_repair: a.id_repair, res_id_original: o, res_id_repair: r, intersect: i, original_except_repair: e1,
                    repair_except_original: e2);

            }).ToList();

            diff = diff.Where(a => a.intersect.Count != a.res_id_original.Count || a.intersect.Count != a.res_id_repair.Count).ToList();

            return diff;
        }

        public static void get_limited_ids()
        {

            var seqs_all = File.ReadAllLines(@"c:\betastrands_dataset\betastrands_dataset_sequences.txt");
            var seqs_limited = File.ReadAllLines(@"c:\betastrands_dataset\betastrands_dataset_sequences_limited.txt");

            var ids = seqs_all.Select((a, i) => seqs_limited.Contains(a) ? i : -1).Distinct().Where(a => a != -1).ToList();

            File.WriteAllLines(@"c:\betastrands_dataset\copy_limited.bat", ids.Select(a=> $@"copy {a}.* limited\").ToList());
            //Console.ReadKey();

        }

        public static bool verbose_debug = false;

        public static void Main(string[] args)
        {
            //var aaa= subsequence_classification_data.aaindex_subset_templates_search();


            //return;
            //r_peptides.init_r();
            //Console.WriteLine(r_peptides.get_values("AAAAAAAAAALLLLLLL"));
            //for (var i = 0; i < 100; i++) Console.ReadLine();

            //return;

            //get_limited_ids();
            //return;

            //var xxxxx = new mpsa_reader(@"c:\betastrands_dataset\ss_meta\0.psipred-swissprot");

            //return;
            ////pssm test
            //var blast_db_folder = @"c:\blast\db\";
            //var blast_dbs =new[]{"nr", "swissprot", "uniref90"};
            //var remotes = new[]{true, false};
            //var run = false;
            //foreach (var blast_db in blast_dbs)
            //{
            //    foreach (var remote in remotes)
            //    {
            //        get_pssms(blast_db, blast_db_folder, remote, run);
            //    }
            //}
            //Console.WriteLine("Finished");
            //Console.ReadLine();
            //return;



            //dhc_dataset_maker.run_dhc_dataset_maker(substructure_type.dimorphic, 1, 2, true, false);
            //return;


            //for (var i = 3; i <= 5; i++)
            //{
            //    var a = new string("ALGYALGY".ToCharArray(),0, i);
            //    var pepValues = r_peptides.get_values(a);
            //    Console.WriteLine($"!!! length {i}: sequence {a}: features: {pepValues.Count}");
            //}

            //Console.ReadKey();
            //return;

            //REngine engine = REngine.GetInstance();
            ////engine.Evaluate(@"install.packages(""Peptides"")");

            //var r_init_cmds = $@"
            //#install.packages(""devtools"", repos=""http://cran.us.r-project.org"")
            //#library(devtools)
            //#install_github(""dosorio/Peptides"")
            //library(Peptides)
            //";
            //r_init_cmds.Split(new char[] { '\r', '\n'}).Where(a=>!string.IsNullOrWhiteSpace(a) && !a.Trim().StartsWith("#")).ToList().ForEach(a=> engine.Evaluate(a));

            //var xyz_e = engine.Evaluate("xyz <- lengthpep(seq = \"GLPRKILCAIAKKKGKCKGPLKLVCKC\")").AsDataFrame();
            //NumericVector xyz = engine.GetSymbol("xyz").AsNumeric();



            //var edge_files = Directory.GetFiles(@"C:\Users\aaron\Desktop\ring_pdb\pdb\", "*.edges");
            ////var node_files = Directory.GetFiles(@"C:\Users\aaron\Desktop\ring_pdb\pdb\","*.nodes");


            //var ring_edges = edge_files.SelectMany(a=> ring.ring_edge.load(a));

            //var i_types0 = ring_edges.SelectMany(a => new[]
            //{
            //    a.Interaction.text
            //}).ToList();

            //var ix = i_types0.Distinct().OrderBy(a => a).Select(a=>(a,i_types0.Count(b=>b==a))).ToList();

            //////var i_types1 = ring_edges.Select(a => a.Interaction.subtype_node1).Distinct().ToList();
            //////var i_types2 = ring_edges.Select(a => a.Interaction.subtype_node2).Distinct().ToList();

            //Console.WriteLine(string.Join("\r\n", ix));

            //SC MC





            //Console.ReadKey();
            //return;
            //compare_repaired_pdb_to_pdb();
            //return;

            //check_ala_scanning_files_are_complete();
            //check_foldx_position_scan_files_are_completed();
            //check_buildmodel_position_scan_files_are_complete();
            //check_buildmodel_subseq_mutant_files_are_complete();
            //Console.ReadLine();
            //return;


            // sable test
            //sable.load(@"c:\sable\test.txt");
            // return;

            // uniprot test
            //get_uniprot_sequences();
            //return;

            // pdb sequences test
            //get_pdb_sequences(limited_to_dimorphics_and_standard: false);
            //get_pdb_sequences(limited_to_dimorphics_and_standard: true);
            //Console.WriteLine("Finished");
            //Console.ReadLine();
            //return;

            var sw1 = new Stopwatch();
            sw1.Start();

            //var fghx = feature_calcs.feature_pse_aac("CCCCCCCCCCCCCEEEEEEEEHHHHHHT", feature_calcs.seq_type.secondary_structure_sequence, new feature_calcs.pse_aac_options(), false, true);
            //Console.WriteLine();
            //get_pdb_sequences();
            //return;
            //var xxx = feature_calcs.feature_pse_aac("AAA", false, false);

            //return;


            //fill_missing_chains();
            //get_missing_data();
            //return;

            //extract_all_pdbs();
            //repair_all_extracted_pdbs();

            //get_pdb_sequences();

            //get_pssms();
            //load_pssms();

            //return;
            // unknown classes
            //var run_larks_unknown = true;
            //if (run_dimorphics_pos && run_dhc_pos) throw new Exception("Is that right?");
            //if (run_coils_neg && run_strands_neg) throw new Exception("Is that right?");
            //var uid = 2;
            //class_list.Add((-1, substructure_type.coil, 3, nameof(substructure_type.coil), $@"c:\bioinf\e_{nameof(substructure_type.coil)}_{uid}", $@"c:\bioinf\d_{nameof(substructure_type.coil)}_{uid}"));
            //class_list.Add((-1, substructure_type.standard_beta_strand, 2, nameof(substructure_type.standard_beta_strand), $@"c:\bioinf\e_{nameof(substructure_type.standard_beta_strand)}_{uid}", $@"c:\bioinf\d_{nameof(substructure_type.standard_beta_strand)}_{uid}"));
            //class_list.Add((-1, substructure_type.standard_beta_strand_with_host_coil, 3, nameof(substructure_type.standard_beta_strand_with_host_coil), $@"c:\bioinf\e_{nameof(substructure_type.standard_beta_strand_with_host_coil)}_{uid}", $@"c:\bioinf\d_{nameof(substructure_type.standard_beta_strand_with_host_coil)}_{uid}"));
            //class_list.Add((-1, substructure_type.standard_beta_strand_full_protein_sequence, 3, nameof(substructure_type.standard_beta_strand_full_protein_sequence), $@"c:\bioinf\e_{nameof(substructure_type.standard_beta_strand_full_protein_sequence)}_{uid}", $@"c:\bioinf\d_{nameof(substructure_type.standard_beta_strand_full_protein_sequence)}_{uid}"));
            //class_list.Add((+1, substructure_type.dimorphics, 2, nameof(substructure_type.dimorphics), $@"c:\bioinf\e_{nameof(substructure_type.dimorphics)}_{uid}", $@"c:\bioinf\d_{nameof(substructure_type.dimorphics)}_{uid}"));
            //class_list.Add((+1, substructure_type.dimorphics_with_host_coil, 3, nameof(substructure_type.dimorphics_with_host_coil), $@"c:\bioinf\e_{nameof(substructure_type.dimorphics_with_host_coil)}_{uid}", $@"c:\bioinf\d_{nameof(substructure_type.dimorphics_with_host_coil)}_{uid}"));
            //class_list.Add((+1, substructure_type.dimorphics_full_protein_sequence, 3, nameof(substructure_type.dimorphics_full_protein_sequence), $@"c:\bioinf\e_{nameof(substructure_type.dimorphics_full_protein_sequence)}_{uid}", $@"c:\bioinf\d_{nameof(substructure_type.dimorphics_full_protein_sequence)}_{uid}"));


            var make_distance_matrices = false;// true;
            var encode_features = true;

            //var sequence_only_features = false;
            //var calculate_1d_features = true;
            //var calculate_3d_features = false;
            //var calculate_1d_nh_features = false;
            //var calculate_3d_nh_features = false;
            var use_dssp3 = true;
            var limit_for_testing = false;
            var max_features = 100;//(int)Math.Round((132 * 2) * 1.1);

            if (limit_for_testing)
            {
                Program.WriteLine($"Warning: {nameof(limit_for_testing)} == {limit_for_testing}");
            }

            // betastrands classifier (standard vs. dimorphic strand)
            var run_dimorphics_pos = false; // note: old classifier for paper
            var run_strands_neg = false; // note: old classifier for paper

            // coils classifier (standard vs. dimorphic coil)
            var run_dhc_pos = false;  // true !!!
            var run_coils_neg = true;

            // Note: DNA Binding: Subsequence only. must be true here - changed to subsequence rather than protein (even though it is protein level in reality)

            // Experiments: This is the first output used for the new experiments 18/04/19

            //subsequence_classification_data.feature_types_1d feature_types_subsequence_1d = new subsequence_classification_data.feature_types_1d()
            //{
            //    pse_aac_sequence_classification_data = true,  // 27,283 features
            //    sequence_geometry_classification_data = true,   // 15 features
            //    sable_classification_data = true, // 7 features
            //    dna_binding_prediction_data = true, // 7 features
            //    intrinsically_unordered_data = true, // 1,921 features

            //    mpsa_classification_data_subsequence = false, // 82,993 features
            //    aa_index_classification_data = false, // ,128,161 features                
            //    blast_pssm_subsequence_classification_data = false, // 2,126,881 features

            //    r_peptides = false // R.NET ERROR - Not working - what happened??
            //};

            // Experiments: Second output

            subsequence_classification_data.feature_types_1d feature_types_subsequence_1d = new subsequence_classification_data.feature_types_1d()
            {
                pse_aac_sequence_classification_data = true,  //1. 27,283 features
                sequence_geometry_classification_data = true,   //2. 15 features
                sable_classification_data = true, //3. 7 features
                dna_binding_prediction_data = true, //4. 7 features
                intrinsically_unordered_data = true, //5. 1,921 features

                mpsa_classification_data_subsequence = true, //6. 82,993 features
                aa_index_classification_data = true, //7. before: ,128,161 features, updated: less
                blast_pssm_subsequence_classification_data = true, //8. 2,126,881 features

                r_peptides = true 
                
                
            };


            //feature_types_subsequence_1d = null;

            subsequence_classification_data.feature_types_1d feature_types_neighbourhood_1d = new subsequence_classification_data.feature_types_1d()
            {
                pse_aac_sequence_classification_data = true,
                sequence_geometry_classification_data = true,

                mpsa_classification_data_subsequence = true,
                intrinsically_unordered_data = true,
                aa_index_classification_data = true,
                sable_classification_data = true,
                dna_binding_prediction_data = true,//must be false - protein level only
                blast_pssm_subsequence_classification_data = true,
                r_peptides = false
            };
            feature_types_neighbourhood_1d = null;

            subsequence_classification_data.feature_types_1d feature_types_protein_1d = new subsequence_classification_data.feature_types_1d()
            {
                pse_aac_sequence_classification_data = true,
                sequence_geometry_classification_data = true,

                mpsa_classification_data_subsequence = true,
                intrinsically_unordered_data = true,
                aa_index_classification_data = true,
                sable_classification_data = true,
                dna_binding_prediction_data = true,
                blast_pssm_subsequence_classification_data = true,
                r_peptides = false
            };
            feature_types_protein_1d = null;


            subsequence_classification_data.feature_types_3d feature_types_subsequence_3d = new subsequence_classification_data.feature_types_3d()
            {
                sasa_classification_data = false,
                intramolecular_classification_data = false,
                aa_aa_distances = false,
                foldx_classification_data = false,
                tortuosity_classification_data = false,
                ring_classification_data = false,
                pse_ssc_dssp_classification_data = false,
            };

            feature_types_subsequence_3d = null;

            subsequence_classification_data.feature_types_3d feature_types_neighbourhood_3d = new subsequence_classification_data.feature_types_3d()
            {
                sasa_classification_data = false,
                intramolecular_classification_data = false,
                aa_aa_distances = false,
                foldx_classification_data = true,
                tortuosity_classification_data = false,
                ring_classification_data = false,
                pse_ssc_dssp_classification_data = false,
            };
            feature_types_neighbourhood_3d = null;

            subsequence_classification_data.feature_types_3d feature_types_protein_3d = new subsequence_classification_data.feature_types_3d()
            {
                sasa_classification_data = false,
                intramolecular_classification_data = false,
                aa_aa_distances = false,
                foldx_classification_data = true,
                tortuosity_classification_data = false,
                ring_classification_data = false,
                pse_ssc_dssp_classification_data = false,
            };

            feature_types_protein_3d = null;



            // disable subsequence / neighbourhood / protein features
            //feature_types_subsequence_1d = null;
            //feature_types_neighbourhood_1d = null;
            //feature_types_protein_1d = null;
            //feature_types_subsequence_3d = null;
            //feature_types_neighbourhood_3d = null;
            //feature_types_protein_3d = null;
            

            var class_info_list = new List<class_info>();

            if (run_strands_neg)
            {
                class_info_list.Add(new class_info()
                {
                    class_id = -1,
                    substructure_type = substructure_type.standard_strand,
                    min_subsequence_length = 2,
                    class_name = nameof(substructure_type.standard_strand),
                });
            }

            if (run_dimorphics_pos)
            {
                class_info_list.Add(new class_info()
                {
                    class_id = +1,
                    substructure_type = substructure_type.dimorphic,
                    min_subsequence_length = 2,
                    class_name = nameof(substructure_type.dimorphic),
                });
            }

            if (run_dhc_pos)
            {
                class_info_list.Add(new class_info()
                {
                    class_id = +1,
                    substructure_type = substructure_type.dimorphic_coil,
                    min_subsequence_length = 3,
                    class_name = nameof(substructure_type.dimorphic_coil),
                });
            }

            if (run_coils_neg)
            {
                class_info_list.Add(new class_info()
                {
                    class_id = -1,
                    substructure_type = substructure_type.standard_coil,
                    min_subsequence_length = 3,
                    class_name = nameof(substructure_type.standard_coil),
                });
            }



            var tasks = new List<Task<List<subsequence_classification_data>>>();

            var max_tasks = 2000;//-1;//Environment.ProcessorCount;

            if (max_tasks < 0)
            {
                max_tasks = 2000;//Environment.ProcessorCount * Math.Abs(max_tasks) * 10;
            }

            var data_lists = new List<List<subsequence_classification_data>>();




            foreach (var class_info in class_info_list)
            {
                if (max_tasks <= 0)
                {
                    var x = find_class_data
                    (
                        class_info,
                        feature_types_subsequence_1d,
                        feature_types_neighbourhood_1d,
                        feature_types_protein_1d,
                        feature_types_subsequence_3d,
                        feature_types_neighbourhood_3d,
                        feature_types_protein_3d,
                        make_distance_matrices,
                        encode_features,
                        max_features,
                        use_dssp3,
                        limit_for_testing,
                        max_tasks
                    );

                    data_lists.Add(x);
                }
                else
                {
                    var task = Task.Run(() =>
                    {
                        return find_class_data(
                                                class_info,
                                                feature_types_subsequence_1d,
                                                feature_types_neighbourhood_1d,
                                                feature_types_protein_1d,
                                                feature_types_subsequence_3d,
                                                feature_types_neighbourhood_3d,
                                                feature_types_protein_3d,
                                                make_distance_matrices,
                                                encode_features,
                                                max_features,
                                                use_dssp3,
                                                limit_for_testing,
                                                max_tasks
                                                );
                    });
                    tasks.Add(task);

                    while (tasks.Count(a => !a.IsCompleted) >= max_tasks)
                    {
                        Task.WaitAny(tasks.ToArray<Task>());
                    }
                }
            }



            Task.WaitAll(tasks.ToArray<Task>());

            if (make_distance_matrices)
            {
                data_lists.AddRange(tasks.Select(a => a.Result).ToList());
                var sequence_lists = data_lists.Select(a => a.Select(b => b.aa_subsequence)).ToList();
                var seq_list_all = sequence_lists.SelectMany(a => a).ToList();

                //make_distance_matrix(string.Join("_", class_list.Select(a => a.class_name).ToList()), seq_list_all, true);
            }



            Program.WriteLine($"{nameof(Main)}: Finished. Press any key to exit...");

            sw1.Stop();
            var ts1 = sw1.Elapsed;

            Console.WriteLine($"Finished: ({ts1.Days}d {ts1.Hours}h {ts1.Minutes}m {ts1.Seconds}s)");

            Console.ReadKey();
            Console.ReadKey();
            Console.ReadKey();
            Console.ReadKey();
        }

        public class instance_meta_data
        {
            public string protein_pdb_id;
            public char protein_chain_id;
            public string protein_dimer_type;

            public int subsequence_class_id;
            public string subsequence_parallelism;
            public string subsequence_symmetry_mode;

            public string subsequence_aa_seq;
            public List<(int residue_index, char i_code, char amino_acid)> subsequence_res_ids;
           
            public string subsequence_dssp_monomer;
            public string subsequence_dssp3_monomer;
            public string subsequence_dssp_multimer;
            public string subsequence_dssp3_multimer;
            public List<(string format, string prediction)> ss_predictions;

            public string protein_aa_seq;
        }

        private static List<subsequence_classification_data> find_class_data
            (
            class_info class_info,
            subsequence_classification_data.feature_types_1d feature_types_subsequence_1d,
            subsequence_classification_data.feature_types_1d feature_types_neighbourhood_1d,
            subsequence_classification_data.feature_types_1d feature_types_protein_1d,
            subsequence_classification_data.feature_types_3d feature_types_subsequence_3d,
            subsequence_classification_data.feature_types_3d feature_types_neighbourhood_3d,
            subsequence_classification_data.feature_types_3d feature_types_protein_3d,
            bool make_distance_matrices,
            bool encode_features,
            int max_features,
            bool use_dssp3,
            bool limit_for_testing,
            int max_tasks
            )
        {
            lock (Program._console_lock)
            {
                Console.WriteLine($"{class_info.class_name}: Started: finding class #{class_info.class_id} {class_info.class_name}");
                Console.WriteLine($"{class_info.class_name}: {DateTime.Now.ToLongDateString()} {DateTime.Now.ToLongTimeString()}");
            }

            var sw1 = new Stopwatch();
            sw1.Start();
            var data_list = dhc_dataset_maker.run_dhc_dataset_maker(class_info.substructure_type, class_info.class_id, class_info.min_subsequence_length, use_dssp3, limit_for_testing, max_tasks);
            sw1.Stop();
            var ts1 = sw1.Elapsed;



            lock (Program._console_lock)
            {
                Console.WriteLine($"{class_info.class_name}: {nameof(data_list)}: {data_list.Count} ({ts1.Days}d {ts1.Hours}h {ts1.Minutes}m {ts1.Seconds}s)");
                Console.WriteLine($"{class_info.class_name}: {DateTime.Now.ToLongDateString()} {DateTime.Now.ToLongTimeString()}");
                Console.WriteLine($"{class_info.class_name}: Finished: finding class #{class_info.class_id} {class_info.class_name}");
                Console.WriteLine();
            }

            //return data_list;

            var output_folder = @"F:\dataset\"; //@"C:\betastrands_dataset\svm_features\";

            if (encode_features)
            {
                lock (Program._console_lock)
                {
                    Console.WriteLine($"{class_info.class_name}: Started: encoding and saving class #{class_info.class_id} {class_info.class_name}");
                    Console.WriteLine($"{class_info.class_name}: {DateTime.Now.ToLongDateString()} {DateTime.Now.ToLongTimeString()}");
                }


                int cl;
                int ct;

                lock (Program._console_lock)
                {
                    Console.WriteLine();
                    cl = Console.CursorLeft + 1;
                    ct = Console.CursorTop;
                    Console.WriteLine($"[{new string('.', data_list.Count)}]");
                    Console.WriteLine();
                }

                var all_output_data = new List<string>();
                var comment_rows = new List<string>();
                var comment_headers = new List<string>() // headers starting with '#' are excluded from copying to classification statistics
                {
                    $"row_index",
                    $"{nameof(instance_meta_data.protein_pdb_id)}",
                    $"{nameof(instance_meta_data.protein_chain_id)}",
                    $"{nameof(instance_meta_data.protein_dimer_type)}",
                    $"{nameof(instance_meta_data.subsequence_class_id)}",
                    $"subsequence_{nameof(class_info.class_name)}",
                    $"{nameof(instance_meta_data.subsequence_parallelism)}",
                    $"{nameof(instance_meta_data.subsequence_symmetry_mode)}",
                    $"{nameof(instance_meta_data.subsequence_aa_seq)}_len",
                    $"{nameof(instance_meta_data.subsequence_aa_seq)}",
                    $"#{nameof(instance_meta_data.protein_aa_seq)}_len",
                    $"#{nameof(instance_meta_data.protein_aa_seq)}",
                    $"{nameof(instance_meta_data.subsequence_res_ids)}",
                    $"{nameof(instance_meta_data.subsequence_dssp_monomer)}",
                    $"{nameof(instance_meta_data.subsequence_dssp_multimer)}",
                    $"{nameof(instance_meta_data.subsequence_dssp3_monomer)}",
                    $"{nameof(instance_meta_data.subsequence_dssp3_multimer)}",
                };

                var sw2 = new Stopwatch();
                sw2.Start();

                var total_featues_check = -1;

                List<subsequence_classification_data.feature_info> header_list = null;
                List<string> header_list_str = null;

                //var header_dupe_check = new List<string>();

                foreach (var row in data_list)
                {
                    var row_index = data_list.IndexOf(row);

                    var row_encoded = subsequence_classification_data.encode_subsequence_classification_data_row(row, max_features, feature_types_subsequence_1d, feature_types_neighbourhood_1d, feature_types_protein_1d, feature_types_subsequence_3d, feature_types_neighbourhood_3d, feature_types_protein_3d);
                        
                        

                    if (header_list == null)
                    {
                        header_list = row_encoded.feature_info.Select(a => new subsequence_classification_data.feature_info(a) { feature_value = 0 }).ToList();


                        //header_list_str = header_list.Select(a => $"{a.source}.{a.@group}.{a.member}.{a.perspective}").ToList();
                        //header_list_str =       header_list.Select(a => $"{a.alphabet}.{a.dimension}.{a.category}.{a.source}.{a.@group}.{a.member}.{a.perspective}").ToList();

                        header_list_str = header_list.Select((a, i) => i.ToString()).ToList();

                        var header_list_str_dupe_check = header_list.Select((a, i) => $"{a.alphabet},{a.dimension},{a.category},{a.source},{a.@group},{a.member},{a.perspective}").ToList();
                        var header_list_str_dupe_check_distinct = header_list_str_dupe_check.Distinct().ToList();

                        if (header_list_str_dupe_check.Count != header_list_str_dupe_check_distinct.Count)
                        {
                            var header_list_str_dupe_check_distinct_count = header_list_str_dupe_check_distinct.Select(a => (item: a, count: header_list_str_dupe_check.Count(b => b == a))).OrderByDescending(a=>a.count).ToList();
                            throw new Exception("duplicate headers found");
                        }

                        //header_dupe_check.Add(header_list_str_c_noindex);

                        var header_list_str_c = header_list.Select((a, i) => $"{i.ToString()},{a.alphabet},{a.dimension},{a.category},{a.source},{a.@group},{a.member},{a.perspective}").ToList();
                        header_list_str_c.Insert(0, "fid,alphabet,dimension,category,source,group,member,perspective");

                        Console.WriteLine("Total features: " + header_list.Count);
                        Console.WriteLine();
                        Console.WriteLine("Total dimension: " + header_list.Select(a => a.dimension).Distinct().Count() + ": " + string.Join(", ", header_list.Select(a => a.dimension).Distinct().Take(20).ToList()) + "...");
                        Console.WriteLine("Total alphabet: " + header_list.Select(a => a.alphabet).Distinct().Count() + ": " + string.Join(", ", header_list.Select(a => a.alphabet).Distinct().Take(20).ToList()) + "...");
                        Console.WriteLine("Total category: " + header_list.Select(a => a.category).Distinct().Count() + ": " + string.Join(", ", header_list.Select(a => a.category).Distinct().Take(20).ToList()) + "...");
                        Console.WriteLine("Total source: " + header_list.Select(a => a.source).Distinct().Count() + ": " + string.Join(", ", header_list.Select(a => a.source).Distinct().Take(20).ToList()) + "...");
                        Console.WriteLine("Total group: " + header_list.Select(a => a.@group).Distinct().Count() + ": " + string.Join(", ", header_list.Select(a => a.@group).Distinct().Take(20).ToList()) + "...");
                        Console.WriteLine("Total member: " + header_list.Select(a => a.member).Distinct().Count() + ": " + string.Join(", ", header_list.Select(a => a.member).Distinct().Take(20).ToList()) + "...");
                        Console.WriteLine("Total perspective: " + header_list.Select(a => a.perspective).Distinct().Count() + ": " + string.Join(", ", header_list.Select(a => a.perspective).Distinct().Take(20).ToList()) + "...");
                        Console.WriteLine();




                        var feature_headers_file = Path.Combine(output_folder, $"h__[{class_info.class_name}].csv");//  $@"c:\betastrands_dataset\svm_features\feature_headers_c{row.class_id}.csv";
                        File.WriteAllLines(feature_headers_file, header_list_str_c);
                        Console.WriteLine("Saved headers: " + feature_headers_file);


                        comment_headers.AddRange(row_encoded.instance_meta_data.ss_predictions.Select(a => a.format).ToList());
                    }

                    var total_features = row_encoded.feature_info.Count;

                    if (total_featues_check > -1 && total_featues_check != total_features)
                    {
                        throw new Exception();
                    }
                    else
                    {
                        total_featues_check = total_features;
                    }

                    var all_feature_values = row_encoded.feature_info.Select(a => a.feature_value.ToString("G17")).ToList();
                    all_output_data.Add(string.Join(",", all_feature_values));



                    var comments = new List<string>();
                    comments.Add(row_index.ToString());
                    comments.Add($"{(row_encoded.instance_meta_data.protein_pdb_id)}");
                    comments.Add($"{(row_encoded.instance_meta_data.protein_chain_id)}");
                    comments.Add($"{(row_encoded.instance_meta_data.protein_dimer_type)}");
                    comments.Add($"{(row_encoded.instance_meta_data.subsequence_class_id)}");
                    comments.Add($"{(class_info.class_name)}");
                    comments.Add($"{(row_encoded.instance_meta_data.subsequence_parallelism)}");
                    comments.Add($"{(row_encoded.instance_meta_data.subsequence_symmetry_mode)}");
                    comments.Add($"{(row_encoded.instance_meta_data.subsequence_aa_seq.Length)}");
                    comments.Add($"{(row_encoded.instance_meta_data.subsequence_aa_seq)}");
                    comments.Add($"#{(row_encoded.instance_meta_data.protein_aa_seq.Length)}");
                    comments.Add($"#{(row_encoded.instance_meta_data.protein_aa_seq)}");
                    comments.Add($"{(string.Join(" ", row_encoded.instance_meta_data.subsequence_res_ids.Select(a => a.amino_acid.ToString() + a.residue_index.ToString() + (a.i_code != default && a.i_code != ' ' ? a.i_code.ToString() : "")).ToList()))}");
                    comments.Add($"{(row_encoded.instance_meta_data.subsequence_dssp_monomer)}");
                    comments.Add($"{(row_encoded.instance_meta_data.subsequence_dssp_multimer)}");
                    comments.Add($"{(row_encoded.instance_meta_data.subsequence_dssp3_monomer)}");
                    comments.Add($"{(row_encoded.instance_meta_data.subsequence_dssp3_multimer)}");



                    comments.AddRange(row_encoded.instance_meta_data.ss_predictions.Select(a => a.prediction).ToList());
                    comment_rows.Add(string.Join(",", comments));

                    lock (Program._console_lock)
                    {
                        var cl2 = Console.CursorLeft;
                        var ct2 = Console.CursorTop;

                        Console.SetCursorPosition(cl, ct);
                        Console.Write("|");
                        ct = Console.CursorTop;
                        cl = Console.CursorLeft;
                        Console.SetCursorPosition(cl2, ct2);
                    }

                }

                sw2.Stop();
                var ts2 = sw2.Elapsed;

                lock (Program._console_lock)
                {
                    Console.WriteLine($"{class_info.class_name}: ({ts2.Days}d {ts2.Hours}h {ts2.Minutes}m {ts2.Seconds}s)");
                    Console.WriteLine($"{class_info.class_name}: {DateTime.Now.ToLongDateString()} {DateTime.Now.ToLongTimeString()}");
                    Console.WriteLine($"{class_info.class_name}: Finished: encoding and saving class #{class_info.class_id} {class_info.class_name}");
                    Console.WriteLine();
                }


                

                // save comments file
                comment_rows.Insert(0, string.Join(",", comment_headers));
                var comment_file = Path.Combine(output_folder, $"c__[{class_info.class_name}].csv");
                Directory.CreateDirectory(Path.GetDirectoryName(comment_file));
                File.WriteAllLines(comment_file, comment_rows);
                Program.WriteLine("Saved: " + comment_file);

                // save features file
                all_output_data.Insert(0, string.Join(",", header_list_str));
                var all_features_file = Path.Combine(output_folder, $"f__[{class_info.class_name}].csv");
                Directory.CreateDirectory(Path.GetDirectoryName(all_features_file));
                File.WriteAllLines(all_features_file, all_output_data);
                Program.WriteLine("Saved: " + all_features_file);

                

            }

            return data_list;
        }

        //public static void make_distance_matrix(string file_marker, List<string> sequence_list, bool cluster = false)
        //{
        //    var aa_matrix = new double[sequence_list.Count, sequence_list.Count];
        //    var aa_group_matrix = new double[sequence_list.Count, sequence_list.Count];
        //    for (var i = 0; i < sequence_list.Count; i++)
        //    {
        //        var seq1 = sequence_list[i];
        //        for (var j = 0; j < sequence_list.Count; j++)
        //        {
        //            if (i < j) continue;

        //            var seq2 = sequence_list[j];

        //            var aa_seq_distance = feature_calcs.aa_sequence_distance(seq1, seq2);
        //            aa_matrix[i, j] = aa_seq_distance;
        //            aa_matrix[j, i] = aa_seq_distance;

        //            var aa_group_seq_distance = feature_calcs.aa_group_sequence_distance(seq1, seq2);
        //            aa_group_matrix[i, j] = aa_group_seq_distance;
        //            aa_group_matrix[j, i] = aa_group_seq_distance;
        //        }
        //    }

        //    if (cluster)
        //    {
        //        Console.WriteLine($"clustering {nameof(aa_matrix)}");
        //        cluster_distance_matrix(aa_matrix);

        //        Console.ReadKey();

        //        Console.WriteLine($"clustering {nameof(aa_group_matrix)}");
        //        cluster_distance_matrix(aa_group_matrix);
        //        Console.ReadKey();

        //    }

        //    File.WriteAllLines($@"c:\bioinf\dm_{file_marker}_{nameof(aa_matrix)}.txt", aa_matrix.ToJagged().Select(a => string.Join(",", a)).ToList());
        //    File.WriteAllLines($@"c:\bioinf\dm_{file_marker}_{nameof(aa_group_matrix)}.txt", aa_group_matrix.ToJagged().Select(a => string.Join(",", a)).ToList());
        //}

        public static List<(int x, int y, double d)> matrix_to_list(double[,] matrix)
        {
            var matrix_list = new List<(int x, int y, double d)>();
            for (var x = 0; x < matrix.GetLength(0); x++)
            {
                for (var y = 0; y < matrix.GetLength(0); y++)
                {
                    if (x < y) continue;

                    matrix_list.Add((x, y, matrix[x, y]));
                }
            }

            matrix_list = matrix_list.OrderBy(a => a.d).ThenBy(a => a.x).ThenBy(a => a.y).ToList();
            return matrix_list;
        }

        public static void cluster_distance_matrix(double[,] matrix)
        {
            var matrix_list = matrix_to_list(matrix);

            var clusters = new List<List<int>>();

            for (var i = 0; i < matrix_list.Count; i++)
            {
                var x_cluster = clusters.First(a => a.Contains(matrix_list[i].x));
                var y_cluster = clusters.First(a => a.Contains(matrix_list[i].y));

                var cluster = new List<int>();
                if (x_cluster != null)
                {
                    cluster.AddRange(x_cluster);
                    clusters.Remove(x_cluster);
                }
                else
                {
                    cluster.Add(matrix_list[i].x);
                }
                if (y_cluster != null)
                {
                    cluster.AddRange(y_cluster);
                    clusters.Remove(y_cluster);
                }
                else
                {
                    cluster.Add(matrix_list[i].y);
                }
                clusters.Add(cluster);

                Console.WriteLine($"{i}: {clusters.Count}: {string.Join(", ", clusters.Select((a, j) => $"{j}: {a.Count}").ToList())}");

            }
        }
    }
}
