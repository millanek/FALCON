from pypeflow.data import PypeLocalFile, makePypeLocalFile, fn
from pypeflow.task import PypeTask, PypeThreadTaskBase, PypeTaskBase
from pypeflow.controller import PypeWorkflow, PypeMPWorkflow, PypeThreadWorkflow
from falcon_kit.FastaReader import FastaReader
from falcon_kit.fc_asm_graph import AsmGraph
import argparse
import glob
import sys
import subprocess as sp
import shlex
import os

def make_dirs(d):
    if not os.path.isdir(d):
        os.makedirs(d)


def get_read_ctg_map(rawread_dir, pread_dir, asm_dir, raw_db, corr_db, asm_sf):

    read_map_dir = os.path.abspath(os.path.join(asm_dir, "read_maps"))
    make_dirs(read_map_dir)

    PypeMPWorkflow.setNumThreadAllowed(12, 12)
    wf = PypeMPWorkflow()

    rawread_db = makePypeLocalFile( os.path.join( rawread_dir, raw_db ) )
    rawread_id_file = makePypeLocalFile( os.path.join( read_map_dir, "raw_read_ids" ) )

    @PypeTask( inputs = {"rawread_db": rawread_db},
               outputs =  {"rawread_id_file": rawread_id_file},
               TaskType = PypeThreadTaskBase,
               URL = "task://localhost/dump_rawread_ids" )
    def dump_rawread_ids(self):
        rawread_db = fn( self.rawread_db )
        rawread_id_file = fn( self.rawread_id_file )
        os.system("DBshow -n %s | tr -d '>' | LD_LIBRARY_PATH= awk '{print $1}' > %s" % (rawread_db, rawread_id_file) )

    wf.addTask( dump_rawread_ids )

    pread_db = makePypeLocalFile( os.path.join( pread_dir, corr_db ) )
    pread_id_file = makePypeLocalFile( os.path.join( read_map_dir, "pread_ids" ) )

    @PypeTask( inputs = {"pread_db": pread_db},
               outputs =  {"pread_id_file": pread_id_file},
               TaskType = PypeThreadTaskBase,
               URL = "task://localhost/dump_pread_ids" )
    def dump_pread_ids(self):
        pread_db = fn( self.pread_db )
        pread_id_file = fn( self.pread_id_file )
        os.system("DBshow -n %s | tr -d '>' | LD_LIBRARY_PATH= awk '{print $1}' > %s" % (pread_db, pread_id_file) )

    wf.addTask( dump_pread_ids )

    wf.refreshTargets() # block

    sg_edges_list = makePypeLocalFile( os.path.join(asm_dir, "sg_edges_list"+asm_sf) )
    utg_data = makePypeLocalFile( os.path.join(asm_dir, "utg_data"+asm_sf) )
    ctg_paths = makePypeLocalFile( os.path.join(asm_dir, "ctg_paths"+asm_sf) )

    inputs = { "rawread_id_file": rawread_id_file,
               "pread_id_file": pread_id_file,
               "sg_edges_list": sg_edges_list,
               "utg_data": utg_data,
               "ctg_paths": ctg_paths }

    read_to_contig_map = makePypeLocalFile( os.path.join(read_map_dir, "read_to_contig_map") )

    @PypeTask( inputs = inputs,
               outputs = {"read_to_contig_map": read_to_contig_map},
               TaskType = PypeThreadTaskBase,
               URL = "task://localhost/get_ctg_read_map" )
    def generate_read_to_ctg_map(self):
        rawread_id_file = fn( self.rawread_id_file )
        pread_id_file = fn( self.pread_id_file )
        read_to_contig_map = fn( self.read_to_contig_map )

        pread_did_to_rid = open(pread_id_file).read().split("\n")
        rid_to_oid = open(rawread_id_file).read().split("\n")

        asm_G = AsmGraph(fn(self.sg_edges_list),
                         fn(self.utg_data),
                         fn(self.ctg_paths) )

        pread_to_contigs = {}

        with open(read_to_contig_map, "w") as f:
            for ctg in asm_G.ctg_data:
                if ctg[-1] == "R":
                    continue
                ctg_g = asm_G.get_sg_for_ctg(ctg)
                for n in ctg_g.nodes():
                    pid = int(n.split(":")[0])

                    rid = pread_did_to_rid[pid].split("/")[1]
                    rid = int(int(rid)/10)
                    oid = rid_to_oid[rid]
                    k = (pid, rid, oid)
                    pread_to_contigs.setdefault( k, set() )
                    pread_to_contigs[ k ].add( ctg )


            for k in pread_to_contigs:
                pid, rid, oid = k
                for ctg in list(pread_to_contigs[ k ]):
                    print >>f, "%09d %09d %s %s" % (pid, rid, oid, ctg)

    wf.addTask( generate_read_to_ctg_map )

    wf.refreshTargets() # block

def parse_args(argv):
    parser = argparse.ArgumentParser(description='generate `2-asm-falcon/read_maps/read_to_contig_map` that contains the \
information from the chain of mapping: (contig id) -> (internal p-read id) -> (internal raw-read id) -> (original read id)',
              formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--basedir', type=str, default="./", help='the base working dir of a FALCON assembly')
    parser.add_argument('--rawReadDir', type=str, default="0-rawreads", help='the (sub)directory where the DAZZLER database with raw PacBio reads is located')
    parser.add_argument('--correctedReadDir', type=str, default="1-preads_ovl", help='the (sub)directory where the DAZZLER database with the corrected PacBio reads is located')
    parser.add_argument('--assemblyDir', type=str, default="2-asm-falcon", help='the (sub)directory with the FALCON assembly contigs')
    parser.add_argument('--rawReadDbName', type=str, default="raw_reads.db", help='the the name of the ACTUAL DAZZLER database with raw PacBio reads')
    parser.add_argument('--correctReadDbName', type=str, default="preads.db", help='the the name of the ACTUAL DAZZLER database with the corrected PacBio reads')
    parser.add_argument('--assemblySuffix', type=str, default="", help='suffix for the assembly filenames sg_edges_list, utg_data, ctg_paths (can specify e.g. the parameters of the assembly)') 
    args = parser.parse_args(argv[1:])
    return args

def main(argv=sys.argv):
    args = parse_args(argv)
    basedir = args.basedir
    rawsubdir = args.rawReadDir
    correctsubdir = args.correctedReadDir
    assembsubdir = args.assemblyDir
    rawread_dir = os.path.abspath( os.path.join( basedir, rawsubdir )  )
    pread_dir = os.path.abspath( os.path.join( basedir, correctsubdir ) )
    asm_dir = os.path.abspath( os.path.join( basedir, assembsubdir) )

    get_read_ctg_map(rawread_dir=rawread_dir, pread_dir=pread_dir, asm_dir=asm_dir, raw_db=args.rawReadDbName, corr_db=args.correctReadDbName, asm_sf=args.assemblySuffix)

if __name__ == "__main__":
    main()
