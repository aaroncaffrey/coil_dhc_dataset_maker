using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Coil_DHC_Dataset_Maker
{
    public class ring
    {

        public class ring_edge
        {
            //http://protein.bio.unipd.it/ring/help
            //"NodeId1     Interaction   NodeId2      Distance   Angle    Energy   Atom1   Atom2   Donor        Positive	Cation	Orientation"
            //A:31:_:SER   HBOND:SC_SC   A:43:_:LYS   2.911      12.333   17.000   OG      NZ      A:43:_:LYS   		

            public (string text, char chain, int res_id, char icode, char amino_acid1, string amino_acid3) NodeId1; // These two columns report the source and the target node. Nodes can be either a residue or a ligand molecule. The node ID format follows the RIN Analyzer and RING (first version) standards.  <chain> : <index> : <insertion_code> : <residue_3_letter_code> 
            public (string text, string interaction_type, string subtype_node1, string subtype_node2) Interaction; // This attribute reports the type of interaction and the interaction subtype. <interaction_type> : <subtype_node1> _ <subtype_node2> Where subtypes values are: main chain (MC), side chain (SC) and ligand (LIG).
            public (string text, char chain, int res_id, char icode, char amino_acid1, string amino_acid3) NodeId2;
            public double Distance;//The distance in Å between atom centers / mass centers / barycenters depending on the type of interaction and the type of residue.
            public double? Angle;//The angle in degree. Depending on the type of interaction it is calculated differently as described in the previous section. For VDW interactions and IAC it is not possible to infer any angle and the -999.9 (NULL) value is assigned.

            public double Energy;//The average bond free energy in Kj/mol according to literature. For hydrogen bonds it varies according to the involved donor/acceptor atoms. For disulphide bonds dissociation enthalpy is reported. For VDW the enrgy corresponds to the average attractiveness component of the van der Waals force. All energy values reported by RING-2.0 are listed in the following table.
            public string Atom1;
            public string Atom2;
            public (string text, char chain, int res_id, char icode, char amino_acid1, string amino_acid3) Donor;
            public string Positive;
            public string Cation;
            public string Orientation;


            public static List<ring_edge> load(string filename)
            {
                var result = new List<ring_edge>();

                var data = File.ReadAllLines(filename).Skip(1).Where(a => !string.IsNullOrWhiteSpace(a)).Select(a=>a.Split('\t').Select(b => b.Trim()).Select(b => b == "-999.9" ? "0" : b).ToList()).ToList();

                foreach (var d in data)
                {

                    var NodeId1_a = (d[0], d[0].Split(':')[0][0], int.Parse(d[0].Split(':')[1]), d[0].Split(':')[2][0], Atom.Aa3To1(d[0].Split(':')[3]), d[0].Split(':')[3]);
                    var Interaction_a = (d[1], d[1].Split(':')[0], d[1].Split(':')[1].Split('_')[0], d[1].Split(':')[1].Split('_')[1]);

                    var NodeId2_a = (d[2], d[2].Split(':')[0][0], int.Parse(d[2].Split(':')[1]), d[2].Split(':')[2][0], Atom.Aa3To1(d[2].Split(':')[3]), d[2].Split(':')[3]);
                    var Donor_a = (d[8], d[8].Length > 0 ? d[8].Split(':')[0][0] : ' ', d[8].Length > 0 ? int.Parse(d[8].Split(':')[1]) : 0, d[8].Length > 0 ? d[8].Split(':')[2][0] : ' ', d[8].Length > 0 ? Atom.Aa3To1(d[8].Split(':')[3]) : ' ', d[8].Length > 0 ? d[8].Split(':')[3] : "");

                    var NodeId2_b = (d[0], d[0].Split(':')[0][0], int.Parse(d[0].Split(':')[1]), d[0].Split(':')[2][0], Atom.Aa3To1(d[0].Split(':')[3]), d[0].Split(':')[3]);
                    var Interaction_b = (d[1], d[1].Split(':')[0], d[1].Split(':')[1].Split('_')[1], d[1].Split(':')[1].Split('_')[0]);

                    var NodeId1_b = (d[2], d[2].Split(':')[0][0], int.Parse(d[2].Split(':')[1]), d[2].Split(':')[2][0], Atom.Aa3To1(d[2].Split(':')[3]), d[2].Split(':')[3]);
                    var Donor_b = (d[8], d[8].Length > 0 ? d[8].Split(':')[0][0] : ' ', d[8].Length > 0 ? int.Parse(d[8].Split(':')[1]) : 0, d[8].Length > 0 ? d[8].Split(':')[2][0] : ' ', d[8].Length > 0 ? Atom.Aa3To1(d[8].Split(':')[3]) : ' ', d[8].Length > 0 ? d[8].Split(':')[3] : "");

                    result.Add(new ring_edge()
                    {
                        NodeId1 = NodeId1_a,
                        Interaction = Interaction_a,
                        NodeId2 = NodeId2_a,
                        Distance = d[3].Length > 0 ? double.Parse(d[3]) : 0,
                        Angle = d[4].Length > 0 && d[4] != "-999.9" ? (double?)double.Parse(d[4]) : null,
                        Energy = d[5].Length > 0 ? double.Parse(d[5]) : 0,
                        Atom1 = d[6],
                        Atom2 = d[7],
                        Donor = Donor_a,
                        Positive = d[9],
                        Cation = d[10],
                        Orientation = d[11],
                    });

                    result.Add(new ring_edge()
                    {
                        NodeId1 = NodeId1_b,
                        Interaction = Interaction_b,
                        NodeId2 = NodeId2_b,
                        Distance = d[3].Length > 0 ? double.Parse(d[3]) : 0,
                        Angle = d[4].Length > 0 && d[4] != "-999.9" ? (double?)double.Parse(d[4]) : null,
                        Energy = d[5].Length > 0 ? double.Parse(d[5]) : 0,
                        Atom1 = d[6],
                        Atom2 = d[7],
                        Donor = Donor_b,
                        Positive = d[9],
                        Cation = d[10],
                        Orientation = d[11],
                    });
                }

                return result;
            }
        }

        public class ring_node
        {
            //http://protein.bio.unipd.it/ring/help
            //"NodeId       Chain   Position   Residue   Dssp    Degree   Bfactor_CA    x         y         z         pdbFileName              Rapdf     Tap"
            // A:31:_:SER   A       31         SER       (tab)   1        0.000         12.051    -15.082   -33.581   1AJYA_Repair.pdb#31.A    -26.022   0.315

            public (string text, char chain, int res_id, char icode, char amino_acid1, string amino_acid3) NodeId; //The node ID.See the corresponding edge attributes for more details about the format.
            public char Chain;//Three columns that reports the residue (or ligand) chain, position (PDB index) and residue 3 letters code.
            public int Position;
            public char Residue1;
            public string Residue3;
            public string Dssp;//The secondary structure calculated with the DSSP algorithm [5] re-implemented in RING-2.0. 
            public double Degree;//The node degree, i.e. the number of directly connected nodes
            public double Bfactor_CA;//The B factor of the C-alpha (when available).
            public double x;//The 3D coordinates of the C-alpha as reported in the original PDB.
            public double y;
            public double z;
            public string pdbFileName;//<pdb_file_name>#<index>.<chain> This is necessary to bind RINanylezer/StructureViz with Chimera visualization of the structure.
            public double Rapdf;//(Optional). The RAPDF energy. It is calculated based on statistical potentials, see [2]. Tosatto, S.C.E., 2005. The victor/FRST function for model quality estimation. J. Comput. Biol. 12, 1316–1327. doi:10.1089/cmb.2005.12.1316
            public double Tap;//(Optional). The TAP energy. It is calculated based on statistical potentials, see [3]. Tosatto, S.C.E., Battistutta, R., 2007. TAP score: torsion angle propensity normalization applied to local protein structure evaluation. BMC Bioinformatics 8, 155. doi:10.1186/1471-2105-8-155

            public static List<ring_node> load(string filename)
            {
                var result = new List<ring_node>();

                var data = File.ReadAllLines(filename).Skip(1).Where(a => !string.IsNullOrWhiteSpace(a)).Select(a => a.Split('\t').Select(b => b.Trim()).Select(b => b == "-999.9" ? "0" : b).ToList()).ToList();

                foreach (var d in data)
                {
                    var NodeId = (d[0], d[0].Split(':')[0][0], int.Parse(d[0].Split(':')[1]), d[0].Split(':')[2][0],
                        Atom.Aa3To1(d[0].Split(':')[3]), d[0].Split(':')[3]);

                    result.Add(new ring_node()
                    {
                        NodeId = NodeId,
                        Chain= d[1][0],
                        Position= d[2].Length > 0 ? int.Parse(d[2]): 0,
                        Residue1 = Atom.Aa3To1(d[3]),
                        Residue3 = d[3],
                        Dssp = d[4],
                        Degree= d[5].Length > 0 ? double.Parse(d[5]): 0,
                        Bfactor_CA = d[6].Length > 0 ? double.Parse(d[6]): 0,
                        x = d[7].Length > 0 ? double.Parse(d[7]): 0,
                        y = d[8].Length > 0 ? double.Parse(d[8]): 0,
                        z = d[9].Length > 0 ? double.Parse(d[9]): 0,
                        pdbFileName = d[10],
                        Rapdf= d[11].Length > 0 ? double.Parse(d[11]): 0,
                        Tap = d[12].Length > 0 ? double.Parse(d[12]): 0,
                    });
                }

                return result;
            }
        }



    }
}
