""" """
import logging
import os
import sys

import numpy as np
import pandas as pd
import torch

from rhofold.data.balstn import BLASTN
from rhofold.rhofold import RhoFold
from rhofold.config import rhofold_config
from rhofold.utils import get_device, save_ss2ct, timing
from rhofold.relax.relax import AmberRelaxation
from rhofold.utils.alphabet import get_features

import tqdm 

def fastawriter(sequence, outputpath):

    fastafile = outputpath+'/seq.fasta'
    with open(fastafile, 'w') as file:
        file.write('> \n')
        file.write(sequence + '\n')
    
    return fastafile

@torch.no_grad()
def main(config):
    '''
    RhoFold Inference pipeline 
    ~ changes occured from IOANNIS.durdagilab for this inference pipeline to work alone with the designed NSGA-II ~
    ~ for original RhoFold+ predictions, please use the inference .py file offered in the official RhoFold+ repo ~
    '''

    logger = logging.getLogger('RhoFold Inference')

    model = RhoFold(rhofold_config)

    model.load_state_dict(torch.load(config.ckpt, map_location=torch.device('cpu'))['model'])
    model.eval()

    config.device = get_device(config.device)
    model = model.to(config.device)

    searchpool = pd.read_csv(config.input_df)

    if config.single_seq_pred and config.input_fas is None:
        
        for index, row in tqdm.tqdm(list(searchpool.iterrows()), desc=f'Predicting structures'):

            idx = row['ID']
            seq = row['Sequence']
            
            # output path
            seqpath = config.output_dir + f'/{idx}' 
            #os.makedirs(seqpath, exist_ok=True)

            fas_path = seqpath+'/aptamer.fasta'
            config.input_a3m = fas_path
            data_dict = get_features(fas_path, config.input_a3m)

            # Forward pass
            outputs = model(tokens=data_dict['tokens'].to(config.device),
                        rna_fm_tokens=data_dict['rna_fm_tokens'].to(config.device),
                        seq=data_dict['seq'],
                        )

            output = outputs[-1]

            # Secondary structure, .ct format
            ss_prob_map = torch.sigmoid(output['ss'][0, 0]).data.cpu().numpy()
            ss_file = f'{seqpath}/ss.ct'
            save_ss2ct(ss_prob_map, data_dict['seq'], ss_file, threshold=0.5)

            # Dist prob map & Secondary structure prob map, .npz format
            #npz_file = f'{seqpath}/results.npz'
            #np.savez_compressed(npz_file,
            #                    dist_n = torch.softmax(output['n'].squeeze(0), dim=0).data.cpu().numpy(),
            #                    dist_p = torch.softmax(output['p'].squeeze(0), dim=0).data.cpu().numpy(),
            #                    dist_c = torch.softmax(output['c4_'].squeeze(0), dim=0).data.cpu().numpy(),
            #                    ss_prob_map = ss_prob_map,
            #                    plddt = output['plddt'][0].data.cpu().numpy(),
            #                    )

            # Save the prediction
            unrelaxed_model = f'{seqpath}/aptamer.pdb'

            # The last cords prediction
            node_cords_pred = output['cord_tns_pred'][-1].squeeze(0)
            model.structure_module.converter.export_pdb_file(data_dict['seq'],
                                                            node_cords_pred.data.cpu().numpy(),
                                                            path=unrelaxed_model, chain_id=None,
                                                            confidence=output['plddt'][0].data.cpu().numpy(),
                                                            logger=logger)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()

    parser.add_argument("--device", help="Default cpu. If GPUs are available, you can set --device cuda:<GPU_index> for faster prediction.", default=None)
    parser.add_argument("--ckpt", help="Path to the pretrained model, default ./pretrained/model_20221010_params.pt", default='./pretrained/rhofold_pretrained_params.pt')
    parser.add_argument("--input_fas", help="Path to the input fasta file. Valid nucleic acids in RNA sequence: A, U, G, C", default=None)
    
    #### IOANNIS.durdagilab #### 
    parser.add_argument("--input_df", help="List of sequences in string format.", required=True)
    #### IOANNIS.durdagilab ####
    
    parser.add_argument("--input_a3m", help="Path to the input msa file. Default None."
                                            "If --input_a3m is not given (set to None), MSA will be generated automatically. ", default=None)
    parser.add_argument("--output_dir", help="Path to the output dir. "
                                             "3D prediction is saved in .pdb format. "
                                             "Distogram prediction is saved in .npz format. "
                                             "Secondary structure prediction is save in .ct format. ", required=True)
    parser.add_argument("--single_seq_pred", help="Default False. If --single_seq_pred is set to True, "
                                                       "the modeling will run using single sequence only (input_fas)", default=False)
    parser.add_argument("--database_dpath", help="Path to the pretrained model, default ./database", default='./database')
    parser.add_argument("--binary_dpath", help="Path to the pretrained model, default ./rhofold/data/bin", default='./rhofold/data/bin')

    args = parser.parse_args()
    main(args)
