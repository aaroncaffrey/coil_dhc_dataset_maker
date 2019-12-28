using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace Coil_DHC_Dataset_Maker
{

    public class coil_dataset_maker
    {
        public static List<subsequence_classification_data> find_coils(string dimer_type, string pdb_id, List<Atom> pdb_atoms, int coils_class_id = -1, int min_coil_length_required = 2, bool use_dssp3 = true)
        {
            pdb_atoms = pdb_atoms.Where(a => string.Equals(a.pdb_id, pdb_id, StringComparison.InvariantCultureIgnoreCase)).ToList();

            //pdb_atoms = pdb_atoms.Where(a => a.Multimer_DSSP3 == 'C').ToList(); // only need info about coils

            var chain_ids = pdb_atoms.Select(a => a.chain_id).Distinct().ToList();

            var tasks = new List<Task<List<subsequence_classification_data>>>();

            var max_tasks = 2000;//Environment.ProcessorCount * 1000;

            foreach (var loop_chain_id in chain_ids)
            {
                var chain_id = loop_chain_id;

                var task = Task.Run(() =>
                {
                    var pdb_chain_coils = new List<subsequence_classification_data>();

                    var pdb_chain_atoms = pdb_atoms.Where(a => a.chain_id == chain_id).ToList();
                    var pdb_chain_master_atoms = Atom.select_amino_acid_master_atoms(pdb_id, pdb_chain_atoms);

                    var coil_subsequence_atoms = new List<Atom>();
                    var coil_subsequence_master_atoms = new List<Atom>();


                    for (var pdb_chain_master_atoms_index = 0; pdb_chain_master_atoms_index < pdb_chain_master_atoms.Count; pdb_chain_master_atoms_index++)
                    {
                        var master_atom = pdb_chain_master_atoms[pdb_chain_master_atoms_index];

                        var this_res_id = master_atom.residue_index;
                        var last_res_id = pdb_chain_master_atoms_index == 0 ? this_res_id - 1 : pdb_chain_master_atoms[pdb_chain_master_atoms_index - 1].residue_index;
                        var next_res_id = pdb_chain_master_atoms_index == pdb_chain_master_atoms.Count - 1 ? this_res_id + 1 : pdb_chain_master_atoms[pdb_chain_master_atoms_index + 1].residue_index;

                        var res_id_consecutive_to_last = this_res_id == last_res_id + 1;
                        var res_id_consecutive_to_next = this_res_id == next_res_id - 1; ;

                        if ((use_dssp3 && master_atom.multimer_dssp3 == 'C') || (!use_dssp3 && master_atom.multimer_dssp == 'C'))
                        {
                            coil_subsequence_master_atoms.Add(master_atom);
                            coil_subsequence_atoms.AddRange(master_atom.amino_acid_atoms);
                        }

                        if (coil_subsequence_master_atoms.Count > 0 && 
                            ( ((use_dssp3 && master_atom.multimer_dssp3 != 'C') || (!use_dssp3 && master_atom.multimer_dssp != 'C')) || !res_id_consecutive_to_next || pdb_chain_master_atoms_index == pdb_chain_master_atoms.Count - 1))
                        {
                            if (coil_subsequence_master_atoms.Count == pdb_chain_master_atoms.Count) break; // whole chain is a coil - not useful

                            if (coil_subsequence_master_atoms.Count >= min_coil_length_required)
                            {
                                var scd = new subsequence_classification_data
                                {
                                    dimer_type = dimer_type,
                                    class_id = coils_class_id,
                                    pdb_id = pdb_id,
                                    chain_id = chain_id,
                                    res_ids = coil_subsequence_master_atoms.Select(a => (a.residue_index, a.i_code, a.amino_acid)).ToList(),
                                    aa_subsequence = string.Join("",
                                            coil_subsequence_master_atoms.Select(a => a.amino_acid).ToList()),
                                    dssp_multimer_subsequence =
                                        string.Join("",
                                            coil_subsequence_master_atoms.Select(a => a.multimer_dssp).ToList()),
                                    dssp_monomer_subsequence =
                                        string.Join("",
                                            coil_subsequence_master_atoms.Select(a => a.monomer_dssp).ToList()),
                                    stride_multimer_subsequence =
                                        string.Join("",
                                            coil_subsequence_master_atoms.Select(a => a.multimer_stride).ToList()),
                                    stride_monomer_subsequence =
                                        string.Join("",
                                            coil_subsequence_master_atoms.Select(a => a.monomer_stride).ToList()),
                                    pdb_chain_atoms = pdb_chain_atoms,
                                    pdb_chain_master_atoms = pdb_chain_master_atoms,
                                    subsequence_atoms = coil_subsequence_atoms,
                                    subsequence_master_atoms = coil_subsequence_master_atoms
                                };
                                var fx = foldx_caller.calc_energy_differences(scd.pdb_id, scd.chain_id, scd.res_ids, false, subsequence_classification_data.protein_data_sources.subsequence_3d);
                                scd.foldx_energy_differences = fx;

                                //subsequence_data.run_foldx_for_mutant_energy_data();
                                //subsequence_data.calculate_classification_data();

                                //var classification_data = subsequence_classification_data.calculate_classification_data(subsequence_data); //pdb_id, chain_id, coil_seq, pdb_chain_atoms, pdb_chain_master_atoms, coil_master_atoms);
                                pdb_chain_coils.Add(scd);
                            }

                            coil_subsequence_master_atoms = new List<Atom>();
                            coil_subsequence_atoms = new List<Atom>();
                        }
                    }

                    return pdb_chain_coils;
                });
                tasks.Add(task);

                while (tasks.Count(a => !a.IsCompleted) >= max_tasks)
                {
                    Task.WaitAny(tasks.ToArray<Task>());
                }
            }

            Task.WaitAll(tasks.ToArray<Task>());

            var pdb_coils = tasks.SelectMany(a => a.Result).ToList();

            pdb_coils = pdb_coils.GroupBy(scd => (scd.pdb_id, scd.aa_subsequence)).Select(group => group.First()).Distinct().ToList();

            return pdb_coils;
        }




        public static List<subsequence_classification_data> run_coil_dataset_maker(int coils_class_id = -1, int min_coil_length_required = 3, bool use_dssp3 = true, bool limit_for_testing = false, int max_tasks = 2000)
        {
            //bool limit_for_testing = false;

            
            // only used to limit the coils data set to the same pdb ids as the dimorphics (for same distribution as found in nature)
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

            dimorphics_data = dimorphics_data.Where(a => a.class_name == "Single").ToList();

            var pdb_id_list = dimorphics_data.Select(a => (dimer_type:a.dimer_type, pdb_id: a.pdb_id)).Distinct().ToList();

            if (limit_for_testing)
            {
                pdb_id_list = pdb_id_list.Take(10).ToList();
            }
            //var console_lock = new object();

            int cl;
            int ct;

            lock (Program._console_lock)
            {
                Console.WriteLine();
                cl = Console.CursorLeft + 1;
                ct = Console.CursorTop;
                Console.WriteLine($"[{new string('.', pdb_id_list.Count)}]");
                Console.WriteLine();
            }

            var tasks = new List<Task<List<subsequence_classification_data>>>();

            if (max_tasks < 0)
            {
                max_tasks = 2000;//Environment.ProcessorCount * Math.Abs(max_tasks) * 10;
            }

            foreach (var loop_pdb_id in pdb_id_list)
            {
                var pdb_id = loop_pdb_id;

                var task = Task.Run(() =>
                {
                    var pdb_atoms = Atom.load_atoms_pdb(pdb_id.pdb_id, true, true, false, true, true, true).Where(a => a.pdb_model_index == 0).SelectMany(a => a.pdb_model_chain_atoms).ToList();


                    var pdb_coils = find_coils(pdb_id.dimer_type, pdb_id.pdb_id, pdb_atoms);

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

                    return pdb_coils;
                });

                tasks.Add(task);
                while (tasks.Count(a => !a.IsCompleted) >= max_tasks)
                {
                    Task.WaitAny(tasks.ToArray<Task>());
                }
            }

            Task.WaitAll(tasks.ToArray<Task>());

            var coils = new List<subsequence_classification_data>(); // sasa value, total aa with contacts, sequence

            tasks.ForEach(a => coils.AddRange(a.Result));


            coils = coils.Where(a => a.aa_subsequence.Length >= min_coil_length_required).ToList();

            Program.WriteLine();
            Program.WriteLine("Total coils: " + coils.Count);
            Program.WriteLine();

            return coils;
        }
    }
}
