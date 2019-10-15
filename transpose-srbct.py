gene_description = list()
ews_count = 0
bl_count = 0
nb_count = 0
rms_count = 0
test_count = 0

def read_srbct_raw_dataset():

    global gene_description
    global ews_count
    global bl_count
    global nb_count
    global rms_count
    global test_count

    linecount = 0

    with open("SRBCT.txt", 'r') as datafile:
        for line in datafile:

            splitted_line_list = line.split("\t")

            if linecount == 1:
                for word in splitted_line_list:
                    if "EWS" in word:
                        ews_count = ews_count + 1
                    elif "BL" in word:
                        bl_count = bl_count + 1
                    elif "NB" in word:
                        nb_count = nb_count + 1
                    elif "RMS" in word:
                        rms_count = rms_count + 1
                    elif "TEST" in word:
                        test_count = test_count + 1
            elif linecount > 1:
                # gene_description.append(splitted_line_list[1])
                print(splitted_line_list[1])

            linecount = linecount + 1

    print(gene_description)






def main():
    
    read_srbct_raw_dataset()

    return

main()