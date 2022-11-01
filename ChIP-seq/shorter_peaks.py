import pandas as pd

PEAKS_FILENAME = '../macs/macs_peaks.fa'
SUMMITS_FILENAME = '../macs/macs_summits.bed'
SHORTER_PEAKS_FILENAME = '../macs/macs_shorter_peaks.fa'

summits = pd.read_csv(SUMMITS_FILENAME, sep="\t", header=None)

# Summit position
summit_position = summits.iloc[:, 2]

count = 0
with open(SHORTER_PEAKS_FILENAME, mode='w') as fa_shorter_peaks:
    with open(PEAKS_FILENAME) as fa_peaks:
        for line in fa_peaks:
            if line.startswith(">"):
                abs_pos = int(line.split(":")[1].split("-")[0])
                fa_shorter_peaks.write(line)
            else:
                summit_rel_pos = summits.iloc[count, 2] - abs_pos
                if (summit_rel_pos-100 < 0):
                    fa_shorter_peaks.write(line[0:summit_rel_pos+100].strip())
                    # fa_shorter_peaks.write('\n')
                    print(line[0:summit_rel_pos+100].strip())
                    count += 1
                elif (summit_rel_pos+100 > len(line)):
                    fa_shorter_peaks.write(line[summit_rel_pos-100:].strip())
                    # fa_shorter_peaks.write('\n')
                    print(line[summit_rel_pos-100:].strip())
                    count += 1
                elif (summit_rel_pos-100 < 0 and summit_rel_pos+100 > len(line)):
                    fa_shorter_peaks.write(line.strip())
                    # fa_shorter_peaks.write('\n')
                    print(line)
                    count += 1
                else:
                    fa_shorter_peaks.write(line[summit_rel_pos-100:summit_rel_pos+100].strip())
                    # fa_shorter_peaks.write('\n')
                    print(line[summit_rel_pos-100:summit_rel_pos+100])
                    count += 1
                
                fa_shorter_peaks.write('\n')
                print(count)