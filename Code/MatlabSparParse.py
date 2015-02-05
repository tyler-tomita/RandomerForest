#!/usr/bin/python

import argparse

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('-f',type=str,dest='FileName',required=True,help='Spardat file to be converted to a Matlab compatible format')
    parser.add_argument('-l',dest='Labels',action='store_true',help='Are class labels included as the first entry in each row?')
    parser.add_argument('-v',dest='Values',action='store_true',help='Are values included after the colon?')
    parser.add_argument('-s',type=int,dest='StartIdx',default=0,help='Is the first column index 0 or 1?')
    args = parser.parse_args()
    return args

def main():
    args = parse_arguments()
    fin = open(args.FileName,'r')
    Lines = fin.readlines()
    fin.close()
    fout1 = open(args.FileName+'.ml','w')
    offset = 1 - args.StartIdx
    if args.Labels:
        fout2 = open(args.FileName + '.Labels','w')
        for LineNum,Line in enumerate(Lines):
            Elements = Line.split()
            fout2.write('{}\n'.format(Elements[0]))
            for Element in Elements[1:]:
                if args.Values:
                    col,value = Element.split(':')
                else:
                    col,value = (Element,1)
                fout1.write('{}\t{}\t{}\n'.format(LineNum+1,int(col)+offset,value))
        fout2.close()
    else:
        for LineNum,Line in enumerate(Lines):
            for Element in Line.split():
                if args.Values:
                    col,value = Element.split(':')
                else:
                    col,value = (Element,1)
                fout1.write('{}\t{}\t{}\n'.format(LineNum+1,int(col)+1,value))
    fout1.close()

if __name__ == '__main__':
    main()
