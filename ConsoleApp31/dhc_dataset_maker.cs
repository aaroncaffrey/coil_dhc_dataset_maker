using System;
using System.Collections.Generic;
using System.Data;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace Coil_DHC_Dataset_Maker
{
    public class dhc_dataset_maker
    {
        public enum substructure_type
        {
                none,

                dimorphic,
                dimorphic_coil,
                dimorphic_full_protein_sequence,

                standard_strand,
                standard_strand_coil,
                standard_strand_full_protein_sequence,

                standard_coil,
        
        }

        public static List<subsequence_classification_data> run_dhc_dataset_maker(substructure_type strand_type, int class_id = 1, int min_subsequence_length = 2, bool use_dssp3 = true, bool limit_for_testing = false, int max_tasks = 2000)//-1)
        {
            //bool limit_for_testing = true;
            //const int min_dhc_length_required = 3;

            //var aaindex_db = aaindex.Load();
            var class_name_in_file = "";
            var include_host_coil = false;
            var full_protein_seq = false;
            

            switch (strand_type)
            {
                case substructure_type.none:
                    return null;

                case substructure_type.standard_coil:
                    var x = coil_dataset_maker.run_coil_dataset_maker(class_id, min_subsequence_length, use_dssp3, limit_for_testing, max_tasks);
                    return x;

                case substructure_type.dimorphic:
                    class_name_in_file = "Single";
                    break;

                case substructure_type.dimorphic_coil:
                    class_name_in_file = "Single";
                    include_host_coil = true;
                    break;

                case substructure_type.standard_strand:
                    class_name_in_file = "Multiple";
                    break;

                case substructure_type.standard_strand_coil:
                    class_name_in_file = "Multiple";
                    include_host_coil = true;
                    break;

                case substructure_type.dimorphic_full_protein_sequence:
                    class_name_in_file = "Single";
                    full_protein_seq = true;
                    break;

                case substructure_type.standard_strand_full_protein_sequence:
                    class_name_in_file = "Multiple";
                    full_protein_seq = true;
                    break;

                default:
                    throw new Exception();
            }

            

            var dimorphics_data_all = File.ReadAllLines(@"C:\betastrands_dataset\csv\distinct dimorphics list.csv")
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

            //var dimorphics_data = File.ReadAllLines(@"C:\betastrands_dataset\csv\distinct dimorphics list.csv").Skip(1).Select(a => a.Split(',')).ToList();
            var dimorphics_data = dimorphics_data_all.Where(a => a.class_name == class_name_in_file).ToList();

            if (limit_for_testing)
            {
                var test_limit = 1;
                dimorphics_data = dimorphics_data.Take(test_limit).ToList();
            }

            //dimorphics_data.ForEach(a=>Atom.run_dssp(a));

            //return null;

            //dimorphics_data = dimorphics_data.Take(10).ToList();
            //dimorphics_data = dimorphics_data.Where(a => a[6].Length >= min_dimorphic_length_required).ToList();
            //var dimorphics_ids = dimorphics_data.Select(a => a.First().ToUpperInvariant()).Distinct().OrderBy(a=>a).ToList();

            //0   1   2   3   4   5   6   7   8   9
            //Pdb Id  Dimer Type  Class Name  Symmetry Mode   Parallelism Chain   Sequence Res Index Pdb Seq Pdb Seq Len

            //var dhc_features = new List<List<(string id, object value)>>();


            int cl;
            int ct;

            lock (Program._console_lock)
            {
                Console.WriteLine();
                cl = Console.CursorLeft + 1;
                ct = Console.CursorTop;
                Console.WriteLine($"[{new string('.', dimorphics_data.Count)}]");
                Console.WriteLine();
            }
            //var console_lock = new object();

            //Console.SetCursorPosition(0,ct+1);

            var tasks = new List<Task<subsequence_classification_data>>();

            if (max_tasks == 0)
            {
                max_tasks = 2000;//Environment.ProcessorCount * 10;
            }
            else if (max_tasks < 0)
            {
                max_tasks = 2000;//Math.Abs(max_tasks) * Environment.ProcessorCount * 10;
            }


            foreach (var di in dimorphics_data)
            {
                //var dimorphics_interface = di;
                var dimorphics_interface = (di.pdb_id, di.dimer_type, di.class_name, di.symmetry_mode, di.parallelism, di.chain_number, di.strand_seq, di.optional_res_index);

                //if (dimorphics_interface[6].Length < min_dhc_length_required) continue;

                var task = Task.Run(() =>
                {

                    var pdb_atoms = Atom.load_atoms_pdb(dimorphics_interface.pdb_id, true, true, false, true, true, true).Where(a => a.pdb_model_index == 0).SelectMany(a => a.pdb_model_chain_atoms).ToList();


                    var pdb_id = dimorphics_interface.pdb_id;
                    var chain_number = dimorphics_interface.chain_number;

                    var dimer_type = dimorphics_interface.dimer_type;
                    var class_name = dimorphics_interface.class_name;
                    var symmetry_mode = dimorphics_interface.symmetry_mode;
                    var parallelism = dimorphics_interface.parallelism;
                    var strand_seq = dimorphics_interface.strand_seq;
                    var optional_res_index = dimorphics_interface.optional_res_index;

                    //Atom.load_dssp(pdb_id, pdb_atoms);
                    //Atom.load_stride(pdb_id, pdb_atoms);
                    //Atom.load_rsa(pdb_id, pdb_atoms);
                    //Atom.load_foldx_ala_scanning(pdb_atoms);

                    //const bool find_intermolecular = false;
                    //const bool find_intramolecular = true;

                    //if (find_intermolecular) Atom.find_intermolecular_contacts(pdb_id, pdb_atoms, (double)5);
                    //if (find_intramolecular) Atom.find_intramolecular_contacts(pdb_id, pdb_atoms, (double)5);

                    var chain_list_default_order = pdb_atoms.First(a => a.chain_order_in_pdb_file != null && a.chain_order_in_pdb_file.Count > 0).chain_order_in_pdb_file;//pdb_atoms.Select(a => a.chain_id).Distinct().ToList();
                    var pdb_master_atoms = Atom.select_amino_acid_master_atoms(pdb_id, pdb_atoms);
                    var sequences = pdb_master_atoms.GroupBy(a => a.chain_id).Select(a => (chain_id: a.Key, sequence: string.Join("", a.Select(b => b.amino_acid).ToList()))).ToList();
                    var chain_id = chain_list_default_order[chain_number].chain_id;

                    pdb_atoms = pdb_atoms.Where(a => a.chain_id == chain_id).ToList();
                    var pdb_chain_atoms = pdb_atoms.Where(a => a.chain_id == chain_id).ToList();
                    var pdb_chain_master_atoms = pdb_master_atoms.Where(a => a.chain_id == chain_id).ToList();

                    var strand_protein = sequences.First(a => a.chain_id == chain_id).sequence;
                    var strand_array_index = strand_protein.IndexOf(strand_seq, StringComparison.InvariantCulture);

                    var res_ids = pdb_chain_master_atoms.Select(a => a.residue_index).ToList();


                    if (!string.IsNullOrWhiteSpace(optional_res_index))
                    {
                        strand_array_index = res_ids.IndexOf(int.Parse(optional_res_index));

                        if (strand_protein.Substring(strand_array_index, strand_seq.Length) != strand_seq)
                        {
                            var indexes = strand_protein.Select((a, i) => i + strand_seq.Length < strand_protein.Length && strand_protein.Substring(i, strand_seq.Length) == strand_seq ? i : -1).Where(a => a > -1).ToList();
                            strand_array_index = indexes.Select(a => (index1: a, index2: Math.Abs(strand_array_index - a))).OrderBy(a => a.index2).First().index1;

                        }

                    }

                    if (strand_protein.Substring(strand_array_index, strand_seq.Length) != strand_seq)
                    {
                        throw new Exception("Not correct position");
                    }

                    var others = dimorphics_data_all.Where(a => /*a != di &&*/ a.pdb_id == dimorphics_interface.pdb_id && a.chain_number == dimorphics_interface.chain_number).ToList();
                    

                    foreach (var o in others)
                    {
                        // here, set interface marker properties for atoms

                        var other_pdb_id = o.pdb_id;
                        var other_chain_number = o.chain_number;
                        var other_dimer_type = o.dimer_type;
                        var other_class_name = o.class_name;
                        var other_symmetry_mode = o.symmetry_mode;
                        var other_parallelism = o.parallelism;
                        var other_strand_seq = o.strand_seq;
                        var other_optional_res_index = o.optional_res_index;

                        var other_chain_list_default_order = pdb_atoms.First(a => a.chain_order_in_pdb_file != null && a.chain_order_in_pdb_file.Count > 0).chain_order_in_pdb_file;//pdb_atoms.Select(a => a.chain_id).Distinct().ToList();
                        var other_pdb_master_atoms = Atom.select_amino_acid_master_atoms(other_pdb_id, pdb_atoms);
                        var other_sequences = other_pdb_master_atoms.GroupBy(a => a.chain_id).Select(a => (chain_id: a.Key, sequence: string.Join("", a.Select(b => b.amino_acid).ToList()))).ToList();
                        var other_chain_id = other_chain_list_default_order[other_chain_number].chain_id;

                        var other_pdb_chain_atoms = pdb_atoms.Where(a => a.chain_id == other_chain_id).ToList();
                        var other_pdb_chain_master_atoms = other_pdb_master_atoms.Where(a => a.chain_id == other_chain_id).ToList();

                        var other_strand_protein = other_sequences.First(a => a.chain_id == other_chain_id).sequence;
                        var other_strand_array_index = other_strand_protein.IndexOf(other_strand_seq, StringComparison.InvariantCulture);

                        var other_res_ids = other_pdb_chain_master_atoms.Select(a => a.residue_index).ToList();

                        if (!string.IsNullOrWhiteSpace(other_optional_res_index))
                        {
                            other_strand_array_index = other_res_ids.IndexOf(int.Parse(other_optional_res_index));

                            if (other_strand_protein.Substring(other_strand_array_index, other_strand_seq.Length) != other_strand_seq)
                            {
                                var other_indexes = other_strand_protein.Select((a, i) => i + other_strand_seq.Length < other_strand_protein.Length && other_strand_protein.Substring(i, other_strand_seq.Length) == other_strand_seq ? i : -1).Where(a => a > -1).ToList();
                                other_strand_array_index = other_indexes.Select(a => (index1: a, index2: Math.Abs(other_strand_array_index - a))).OrderBy(a => a.index2).First().index1;
                            }

                        }

                        if (other_strand_protein.Substring(other_strand_array_index, other_strand_seq.Length) != other_strand_seq)
                        {
                            throw new Exception("Not correct position");
                        }

                        var dimorphic_dssp3 = "";

                        for (var i = other_strand_array_index; i < other_strand_array_index + other_strand_seq.Length; i++)
                        {
                            other_pdb_chain_master_atoms[i].parallelism = other_parallelism.Trim().First();
                            other_pdb_chain_master_atoms[i].symmetry_mode = other_symmetry_mode.Trim().First();
                            other_pdb_chain_master_atoms[i].strand_interface_type = other_class_name.Trim().First();

                            other_pdb_chain_master_atoms[i].is_strand_interface_atom = 1;
                            other_pdb_chain_master_atoms[i].amino_acid_atoms.ForEach(a => a.is_strand_interface_atom = 1);

                            if (string.Equals(other_class_name, "Single", StringComparison.InvariantCultureIgnoreCase)|| string.Equals(other_class_name, "Dimorphic", StringComparison.InvariantCultureIgnoreCase))
                            {
                                other_pdb_chain_master_atoms[i].is_dimorphic_strand_interface_atom = 1;
                                other_pdb_chain_master_atoms[i].amino_acid_atoms.ForEach(a => a.is_dimorphic_strand_interface_atom = 1);
                                dimorphic_dssp3 += other_pdb_chain_master_atoms[i].monomer_dssp3;
                            }
                            else if (string.Equals(other_class_name, "Multiple", StringComparison.InvariantCultureIgnoreCase)|| string.Equals(other_class_name, "Standard", StringComparison.InvariantCultureIgnoreCase))
                            {
                                other_pdb_chain_master_atoms[i].is_standard_strand_interface_atom = 1;
                                other_pdb_chain_master_atoms[i].amino_acid_atoms.ForEach(a => a.is_standard_strand_interface_atom = 1);
                            }
                            else if (string.Equals(other_class_name, "Hybrid", StringComparison.InvariantCultureIgnoreCase))
                            {
                                other_pdb_chain_master_atoms[i].is_hybrid_strand_interface_atom = 1;
                                other_pdb_chain_master_atoms[i].amino_acid_atoms.ForEach(a => a.is_hybrid_strand_interface_atom = 1);
                            }

                        }
                        //if (!string.IsNullOrWhiteSpace(dimorphic_dssp3)) Console.WriteLine(dimorphic_dssp3);
                    }

  

                    foreach (var g in pdb_atoms.GroupBy(a => (a.pdb_id , a.chain_id)))
                    {
                        var p = g.Key.pdb_id;
                        var c = g.Key.chain_id;
                        var d = g.ToList();

                        Atom.compare_dssp_to_mpsa(p,c,d);
                    }

                

                    //for (var i = strand_array_index; i < strand_array_index + strand_seq.Length; i++)
                    //{
                    //    pdb_chain_master_atoms[i].is_strand_interface_atom = 1;
                    //    pdb_chain_master_atoms[i].amino_acid_atoms.ForEach(a=>a.is_strand_interface_atom=1);
                    //}

                    var strand_sequence_first_index = strand_array_index;
                    var strand_sequence_last_index = (strand_array_index + strand_seq.Length) - 1;
                    //var dimorphic_sequence_centre = (double)(dimorphic_sequence_first + dimorphic_sequence_last) / (double)2;
                    //var dimorphic_sequence_length = dimorphic_seq.Length;

                    var shc_sequence_before = "";
                    var shc_start_array_index = strand_sequence_first_index;

                    if (include_host_coil)
                    {

                        var continuous_strand1 = (use_dssp3 ? pdb_chain_master_atoms[strand_sequence_first_index].multimer_dssp3 : pdb_chain_master_atoms[strand_sequence_first_index].multimer_dssp) == 'E';

                        for (var i = strand_sequence_first_index - 1; i >= 0; i--)
                        {
                            var index_dssp_code = use_dssp3 ? pdb_chain_master_atoms[i].multimer_dssp3 : pdb_chain_master_atoms[i].multimer_dssp;

                            continuous_strand1 = continuous_strand1 && index_dssp_code == 'E';

                            if (!continuous_strand1 && index_dssp_code != 'C') break;

                            shc_start_array_index = i;
                            shc_sequence_before = pdb_chain_master_atoms[i].amino_acid + shc_sequence_before;
                        }
                    }
                    //if (dhc_start_array_index == -1) dhc_start_array_index = 0;


                    var shc_sequence_after = "";

                    var shc_end_array_index = strand_sequence_last_index;

                    if (include_host_coil)
                    {
                        var continuous_strand2 = (use_dssp3 ? pdb_chain_master_atoms[strand_sequence_last_index].multimer_dssp3 : pdb_chain_master_atoms[strand_sequence_last_index].multimer_dssp) == 'E';

                        for (var i = strand_sequence_last_index + 1; i < pdb_chain_master_atoms.Count; i++)
                        {
                            var index_dssp_code = use_dssp3 ? pdb_chain_master_atoms[i].multimer_dssp3 : pdb_chain_master_atoms[i].multimer_dssp;

                            continuous_strand2 = continuous_strand2 && index_dssp_code == 'E';

                            if (!continuous_strand2 && index_dssp_code != 'C') break;

                            shc_end_array_index = i;
                            shc_sequence_after = shc_sequence_after + pdb_chain_master_atoms[i].amino_acid;
                        }
                    }


                    //if (dhc_end_array_index == -1) dhc_end_array_index = pdb_chain_master_atoms.Count() - 1;

                    var shc_subsequence_master_atoms = pdb_chain_master_atoms.Skip(shc_start_array_index).Take((shc_end_array_index - shc_start_array_index) + 1).ToList();

                    if (full_protein_seq)
                    {
                        shc_subsequence_master_atoms = pdb_chain_master_atoms;
                        shc_start_array_index = 0;
                        shc_end_array_index = pdb_chain_master_atoms.Count - 1;
                        strand_seq = string.Join("", shc_subsequence_master_atoms.Select(a => a.amino_acid).ToList());
                        shc_sequence_before = "";
                        shc_sequence_after = "";
                    }

                    var shc_subsequence_atoms = pdb_chain_atoms.Where(a => shc_subsequence_master_atoms.Any(b => a.pdb_id == b.pdb_id && a.chain_id == b.chain_id && a.residue_index == b.residue_index && a.i_code == b.i_code)).ToList();

                    var shc_sequence = string.Join("", shc_subsequence_master_atoms.Select(a => a.amino_acid).ToList());

                    var shc_dssp_multimer_sequence = string.Join("", shc_subsequence_master_atoms.Select(a => a.multimer_dssp).ToList());
                    var shc_dssp_monomer_sequence = string.Join("", shc_subsequence_master_atoms.Select(a => a.monomer_dssp).ToList());

                    var shc_stride_multimer_sequence = string.Join("", shc_subsequence_master_atoms.Select(a => a.multimer_stride).ToList());
                    var shc_stride_monomer_sequence = string.Join("", shc_subsequence_master_atoms.Select(a => a.monomer_stride).ToList());

                    var shc_length = (shc_end_array_index - shc_start_array_index) + 1;

                    if (shc_sequence.Length != shc_length) throw new Exception("shc length error");
                    if (shc_sequence_before + strand_seq + shc_sequence_after != shc_sequence) throw new Exception("shc seq error");

                    //if (shc_sequence.Length < min_strand_host_coil_length) continue;



                    var scd = new subsequence_classification_data
                    {
                        dimer_type = dimer_type,
                        parallelism = parallelism,
                        symmetry_mode = symmetry_mode,
                        class_id = class_id,
                        //class_name = class_name,
                        pdb_id = pdb_id,
                        chain_id = chain_id,
                        res_ids = shc_subsequence_master_atoms.Select(a => (a.residue_index, a.i_code, a.amino_acid)).ToList(),

                        aa_subsequence = shc_sequence,

                        dssp_multimer_subsequence = shc_dssp_multimer_sequence,
                        dssp_monomer_subsequence = shc_dssp_monomer_sequence,

                        stride_multimer_subsequence = shc_stride_multimer_sequence,
                        stride_monomer_subsequence = shc_stride_monomer_sequence,

                        pdb_chain_atoms = pdb_chain_atoms,
                        pdb_chain_master_atoms = pdb_chain_master_atoms,
                        subsequence_atoms = shc_subsequence_atoms,
                        subsequence_master_atoms = shc_subsequence_master_atoms,
                        
                        
                    };

                    var fx = foldx_caller.calc_energy_differences(scd.pdb_id, scd.chain_id, scd.res_ids, false, subsequence_classification_data.protein_data_sources.subsequence_3d);
                    scd.foldx_energy_differences = fx;

                    //subsequence_classification_data.calculate_classification_data(subsequence_data, subsequence_data.subsequence_master_atoms, source, max_features);



                    lock (Program._console_lock)
                    {
                        //Console.WriteLine($"{dimorphic_seq},{subsequence_data.aa_subsequence},{subsequence_data.dssp_multimer_subsequence},{subsequence_data.dssp_monomer_subsequence},{dhc_stride_multimer_sequence},{dhc_stride_monomer_sequence}");
                    
                        var cl2 = Console.CursorLeft;
                        var ct2 = Console.CursorTop;

                        Console.SetCursorPosition(cl, ct);
                        Console.Write("|");
                        ct = Console.CursorTop;
                        cl = Console.CursorLeft;
                        Console.SetCursorPosition(cl2, ct2);
                    }

                    return scd;
                });

                tasks.Add(task);

                while (tasks.Count(a => !a.IsCompleted) >= max_tasks)
                {
                    Task.WaitAny(tasks.ToArray<Task>());
                }
            }


            Task.WaitAll(tasks.ToArray<Task>());
            var subsequence_classification_data_list = new List<subsequence_classification_data>();

            tasks.ForEach(a => subsequence_classification_data_list.Add(a.Result));

            // check min sequence length requirement is met
            subsequence_classification_data_list = subsequence_classification_data_list.Where(a => a.aa_subsequence.Length >= min_subsequence_length).ToList();

            Program.WriteLine();
            return subsequence_classification_data_list;
        }
    }
}