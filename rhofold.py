import subprocess
import warnings
warnings.simplefilter("ignore") 

def Rhofolder(inputdf, outputfolder):
    command = [
        'python', './strpred/inference4nsga2.py',
        '--input_df', inputdf,
        '--single_seq_pred', 'True',
        '--output_dir', outputfolder,
        '--ckpt', './strpred/pretrained/RhoFold_pretrained.pt'
    ]

    try:
        subprocess.run(command, check=True)

    except subprocess.CalledProcessError as e:
        print(f"Error executing command: {e}")

