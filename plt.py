"a script can plot output file to a '.eps' based gnuplot"
import os
import subprocess
import re
import getopt

def plot_with_gnuplot(file):
    a_res,b_res = file[:-4].split("&")
    standard_plot =['set encoding iso_8859_1\n',
                    'set termoption dash\n',
                    'set parametric\n',
                    'set term post eps color solid enh\n',
                    f'set output "{file[:-4]}.eps"\n',
                    'unset key\n',
                    'set xlabel font "Times-Roman,20"\n',
                    'set ylabel font "Times-Roman,20"\n',
                    'set xtics font "Times-Roman,20"\n',
                    'set ytics font "Times-Roman,20"\n',
                    'set key spacing 1\n',
                    'set key font "Times-Roman,14" \n',
                    'set ytics nomirror\n',
                    'set border lw 2\n',
                    'set xlabel offset 0,-0.5\n',
                    'set ylabel offset -1.2, 0\n',
                    'set size ratio 0.3\n',
                    'set key bottom right\n',
                    'set xrange [0:500]\n',
                    'set yrange [0:10]\n',
                    'set xtics 0,100,500\n',
                    'set ytics 0,2, 10\n']
    '''
    for pair in res_pair:
        a_res,b_res = res_num_dict[pair[0]],res_num_dict[pair[1]]
        with open(f"RRCS-{protain}-{a_res}&{b_res}-new.txt",'a') as f:
            for t in range(1,md_time+1):
                f.write('10.4%f10.5%f' % (t,new_dict[a_res][b_res][t] ))
        with open(f"RRCS-{protain}-{a_res}&{b_res}-old.txt",'a') as f:
            for t in range(1,md_time+1):
                f.write('10.4%f10.5%f' % (t,contact_score[a_res][b_res][t] ))
    '''
    with open(f"{a_res}&{b_res}new.plt",'w') as f:

        for line in standard_plot:
            f.write(line)
        standard2_plot = [ f'plot "{file}"    using ($1/1000):($2*10) notitle w l lw 3 lc rgb "#CCE5FF" ,\\\n',
                            f'"{file}"  using ($1/1000):($2*10) smooth bezier  title \'RRCS-{a_res}&{b_res}\' w l lw 4 lt 2 lc rgb "#0000FF"']
        for line in standard2_plot:
            f.write(line)
    gnuplot= subprocess.Popen('gnuplot',stdin = subprocess.PIPE, stderr=subprocess.PIPE,stdout=subprocess.PIPE)
    gnuplot.stdin .write(f"{a_res}&{b_res}new.plt \n".encode('utf-8'))
    gnuplot.stdin.flush()
    gnuplot_cmd = ['gnuplot', '-p', f"{a_res}&{b_res}new.plt"]
    subprocess.run(gnuplot_cmd)

files_list = os.listdir('.')
print(files_list)
for file in files_list:
    if ".xvg" in file:
        print(file)
        plot_with_gnuplot(file)
