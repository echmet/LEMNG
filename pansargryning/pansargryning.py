#! /usr/bin/env python3

import argparse
import enum
import json
import os
import subprocess
import sys


def _list_to_str(array, formatter, last_formatter=None):
    s = ''

    if array:
        for idx in range(0, len(array) - 1):
            s += formatter(array[idx])
        if last_formatter is None:
            s += formatter(array[-1])
        else:
            s += last_formatter(array[-1])

    return s


def _list_walk(l, handler, handler_last=None):
    if l:
        for idx in range(0, len(l) - 1):
            handler(l[idx])
        if handler_last is None:
            handler(l[-1])
        else:
            handler_last(l[-1])

def _list_to_initer_list(l):
    def fmtr(arg):
        return '{}, '.format(arg)

    def last_fmtr(arg):
        return '{}'.format(arg)

    return '{ ' + _list_to_str(l, fmtr, last_fmtr) + ' }'

class CType(enum.Enum):
    NUCLEUS = 1,
    LIGAND = 2


class CProgram:
    def __init__(self):
        self._includes = []
        self._prelude = []
        self._functions = []

    def add_include(self, code):
        self._includes.append(code)

    def add_prelude(self, code):
        self._prelude.append(code)

    def add_function(self, f):
        self._functions.append(f)

    def __str__(self):
        s = ''

        for i in self._includes:
            s += '#include {}\n'.format(i)

        s += '\n\n'

        for i in self._prelude:
            s += i + '\n'

        s += '\n\n'

        for f in self._functions:
            s += '{}\n\n'.format(f)

        return s


class CFunctionArg:
    def __init__(self, argtype, argname):
        self.argtype = argtype
        self.argname = argname


class CFunction:
    def __init__(self, rtype, name, args):
        self._rtype = rtype
        self._name = name
        self._args = args
        self._code_blocks = []

    def __str__(self):
        def mk_args():
            s = ''

            def fmtr(arg):
                return '{} {}, '.format(arg.argtype, arg.argname)

            def last_fmtr(arg):
                return '{} {}'.format(arg.argtype, arg.argname)

            s += _list_to_str(self._args, fmtr, last_fmtr)

            return s

        s = '{} {}({})\n{{\n'.format(self._rtype, self._name, mk_args())

        for block in self._code_blocks:
            s += '{}\n\n'.format(block)

        s += '}'

        return s

    def add_block(self, block):
        self._code_blocks.append(block.__str__(1))


class CBlock:
    _SPACER = '   '

    @staticmethod
    def make(items):
        block = CBlock()
        for i in items:
            block.add_line(i)
        return block

    def __init__(self):
        self._elems = []

# TODO: Rename this to something sane!!!
    def add_line(self, elem):
        self._elems.append(elem)

    @staticmethod
    def _mk_tabs(indent):
        spcs = ''
        for _ in range(0, indent):
            spcs += CBlock._SPACER

        return spcs

    def __str__(self, indent=0):
        spcs = self._mk_tabs(indent)

        def strify(obj):
            if isinstance(obj, CBlock):
                return obj.__str__(indent+1)
            return '{}{}'.format(spcs, obj)

        def fmtr(item):
            return strify(item) + '\n'

        def fmtr_last(item):
            return strify(item)

        return _list_to_str(self._elems, fmtr, fmtr_last)


class StructInitializer(CBlock):
    def __init__(self, stype, name):
        self._stype = stype
        self._name = name
        self._items = []

    def __str__(self, indent=0):
        tabs = self._mk_tabs(indent)

        s = '{}{} {}{{\n'.format(tabs, self._stype, self._name)

        def fmtr(code):
            return '\t{}{},\n'.format(tabs, code)

        def last_fmtr(code):
            return '\t{}{}\n'.format(tabs, code)

        s += _list_to_str(self._items, fmtr, last_fmtr)

        s += '{}}};'.format(tabs)

        return s

    def add_item(self, code):
        self._items.append(code)


def cname_from_name(name):
    return name.lower().replace(' ', '_')


def gen_header():
    s = '#include <cstdlib>\n'
    s += '#include <iostream>\n\n'
    s += '#include \"barsarkagang_tests.h\"\n\n'
    s += 'using namespace ECHMET;\n'
    s += 'using namespace ECHMET::Barsarkagang;\n'


def gen_ctype(ctype):
    if ctype == CType.NUCLEUS:
        return 'SysComp::ConstituentType::NUCLEUS'
    elif ctype == CType.LIGAND:
        return 'SysComp::ConstituentType::LIGAND'
    raise Exception("Invalid constituent type")


def gen_fixedString(s):
    return 'createFixedString(\"{}\")'.format(s)


def gen_realvec(items):
    s = 'mkRealVec({ '
    s += _list_to_initer_list(items)
    s += ' })'

    return s


def gen_complexes(cpxs, name):
    if cpxs is None:
        return 'nullptr'

    if not cpxs:
        return 'noComplexes()'

    genfunc = CFunction('InCFVec', 'gen_complexforms_' + cname_from_name(name), [])
    outerblock = CBlock()
    outerblock.add_line('{')

    for incf in cpxs:
        ocfblock = CBlock.make(['/* InCFVec */', '{'])

        cfblock = CBlock()
        cfblock.add_line('{},'.format(incf['nucleusCharge']))

        cfblock.add_line('/* InLGVec */')
        for inlg in incf['ligandGroups']:
            cfblock.add_line('{')

            lgblock = CBlock()

            lgblock.add_line('/* InLigandForm */')
            for inlf in inlg['ligands']:
                lgblock.add_line('{')

                lfblock = CBlock()
                lfblock.add_line('\"{}\",'.format(inlf['name']))
                lfblock.add_line('{},'.format(inlf['charge']))
                lfblock.add_line('{},'.format(inlf['maxCount']))
                lfblock.add_line('{},'.format(_list_to_initer_list(inlf['pBs'])))
                lfblock.add_line('{}'.format(_list_to_initer_list(inlf['mobilities'])))

                lgblock.add_line(lfblock)
                lgblock.add_line('}')

            cfblock.add_line(lgblock)
            cfblock.add_line('},')

        ocfblock.add_line(cfblock)
        ocfblock.add_line('},')

        outerblock.add_line(ocfblock)

    outerblock.add_line('}')

    genfunc.add_block(outerblock)

    print(genfunc)

    return 'XXX'


def gen_constituent(name, ctype, chargeLow, chargeHigh, pKas, mobilities, viscosity,
                    complexations):
    ctuent = StructInitializer('SysComp::InConstituent', cname_from_name(name))

    ctuent.add_item(gen_ctype(ctype))
    ctuent.add_item(gen_fixedString(name))
    ctuent.add_item(chargeLow)
    ctuent.add_item(chargeHigh)
    ctuent.add_item(gen_realvec(pKas))
    ctuent.add_item(gen_realvec(mobilities))
    ctuent.add_item(complexations)
    ctuent.add_item(viscosity)

    return ctuent


def gen_complexMap(name, m):
    block = CBlock.make(['CMapping {} = {{'.format(name)])

    e = enumerate(m)

    for idx, k in e:
        v = m[k]
        if idx == len(m) - 1:
            block.add_line('\t{{ \"{}\", {} }}'.format(k, v))
        else:
            block.add_line('\t{{ \"{}\", {} }},'.format(k, v))

    block.add_line('};')

    return block


def gen_calculate(BGElist, sampleList):
    block = CBlock.make(['const auto r = calculate('])

    def mk_vec(l):
        block.add_line('\t{')

        def hndl(item):
            block.add_line('\t\t{},'.format(cname_from_name(item)))

        def hndl_last(item):
            block.add_line('\t\t{}'.format(cname_from_name(item)))

        _list_walk(l, hndl, hndl_last)

        block.add_line('\t},')

    mk_vec(BGElist)
    mk_vec(sampleList)

    block.add_line('\tcBGE, cSample);')

    return block


def gen_check_BGE(expected):
    return CBlock.make(['checkBGE(r, {}, {}, {}, {});'.format(
        expected[0], expected[1], expected[2], expected[3])])


def gen_check_eigenzone(expected):
    return CBlock.make(['checkEigenzone(r.eigenzone, {}, {}, {}, {});'.format(
        expected[0], expected[1], expected[2], expected[3])])

def process_input(root):
    ctuents = root['constituents']

    c_ctuents = []
    concsBGE = dict()
    concsSample = dict()
    BGElist = []
    sampleList = []

    for c in ctuents:
        def mk_cpx():
            if 'complexForms' in c:
                return gen_complexes(c['complexForms'], c['name'])
            return None

        def mk_ctype():
            if c['type'] == 'N':
                return CType.NUCLEUS
            elif c['type'] == 'L':
                return CType.LIGAND
            raise Exception('Invalid constituent type in JSON')

        cct = gen_constituent(c['name'],
                              mk_ctype(),
                              c['chargeLow'],
                              c['chargeHigh'],
                              c['pKas'],
                              c['mobilities'],
                              c['viscosityCoefficient'],
                              mk_cpx())
        c_ctuents.append(cct)

        concsBGE[c['name']] = c['concentrationBGE']
        concsSample[c['name']] = c['concentrationSample']
        sampleList.append(c['name'])
        if c['role'] == 'B':
            BGElist.append(c['name'])

    return (c_ctuents, concsBGE, concsSample, BGElist, sampleList)


def results_file(infile):
    los = infile.rfind('/')
    lod = infile.rfind('.')

    return infile[los+1:lod] + '_results.txt'


def get_expected_results(genpath, infile, ecl_path, lemng_path, is_corr):
    genpath_abs = os.path.abspath(genpath)
    infile_abs = os.path.abspath(infile)
    ecl_path_abs = os.path.abspath(ecl_path)
    lemng_path_abs = os.path.abspath(lemng_path)

    def read_BGE(lines):
        """ Read BGE properties. They are expected to be in the following order
            1) pH 2) Conductivity 3) Ionic strength 4) Buffer capacity
        """
        data = []

        for _ in range(0, 4):
            line = lines.pop(0)
            data.append(float(line))
        line = lines.pop(0)
        if line != '':
            raise Exception('Invalid ouput from results generator - BGE separator')

        if len(data) != 4:
            raise Exception('Invalid ouput from results generator - BGE: insufficient data')

        return data

    def read_eigenzones(lines):
        """ Read eigenzone properties. Each eigenzone is supposed to be separated by an empty line.
            Order of items in an eigenzone block is as follows:
            1) mobility 2) uEMD 3) pH 4) Conductivity
        """

        data = []
        while len(lines) > 0:
            ez = []
            for _ in range(0, 4):
                line = lines.pop(0)
                ez.append(float(line))
            line = lines.pop(0)
            if line != '':
                raise Exception('Invalid ouput from results generator - Eigenzones')

            data.append(ez)

        return data

    resfile = os.getcwd() + '/' + results_file(infile_abs)

    current_lpath = '' if 'LD_LIBRARY_PATH' not in os.environ else os.environ["LD_LIBRARY_PATH"]

    os.environ['LD_LIBRARY_PATH'] = '$LD_LIBRARY_PATH:{}:{}'.format(ecl_path_abs,
                                                                    lemng_path_abs)

    params = [genpath_abs, infile_abs, resfile]

    def extend_switch(p, b):
        p.extend([str(int(b))])
        return p

    params = extend_switch(params, is_corr & 1)
    params = extend_switch(params, is_corr & 2)
    params = extend_switch(params, is_corr & 4)

    try:
        ret = subprocess.run(params)
        if ret.returncode != 0:
            raise Exception('Failed to calculate reference data')
    except Exception as ex:
        raise ex
    finally:
        os.environ['LD_LIBRARY_PATH'] = current_lpath

    data = []
    fh = open(resfile, 'r')

    lines = fh.readlines()
    os.unlink(resfile)

    lines = [l.replace('\n', '') for l in lines]

    data.append(read_BGE(lines))
    for ez in read_eigenzones(lines):
        data.append(ez)

    return data


def make_argparser():
    parser = argparse.ArgumentParser(description='LEMNG Unit tests generator')
    parser.add_argument('--input', help='Input JSON file', type=str)
    parser.add_argument('--output', help='Output C++ file', type=str)
    parser.add_argument('--generator_path', help='Path to reference results generator', type=str)
    parser.add_argument('--ECL_path', help='Path to ECHMETCoreLibs libraries binaries', type=str)
    parser.add_argument('--LEMNG_path', help='Path to LEMNG library binary', type=str)
    parser.add_argument('--debhue', help='Enable Debye-Hückel correction', action='store_true')
    parser.add_argument('--onsfuo', help='Enable Onsager-Fuoss correction', action='store_true')
    parser.add_argument('--viscos', help='Enable viscosity correction', action='store_true')

    parser.set_defaults(debhue=False)
    parser.set_defaults(onsfuo=False)
    parser.set_defaults(viscos=False)

    return parser


def main(outfile, infile, genpath, ecl_path, lemng_path, debhue, onsfuo,
         viscos):
    is_corr = 0
    if debhue:
        is_corr += 1
    if onsfuo:
        is_corr += 2
    if viscos:
        is_corr += 4

    expected = get_expected_results(genpath, infile, ecl_path, lemng_path,
                                    is_corr)

    fh = open(infile, 'r')
    sysdef = json.load(fh)

    c_list, concsBGE, concsSample, BGElist, sampleList = process_input(sysdef)

    cmain = CFunction('int', 'main', [CFunctionArg('int', ''),
                                      CFunctionArg('char **', '')])

    for c in c_list:
        cmain.add_block(c)

    cmain.add_block(gen_complexMap('cBGE', concsBGE))
    cmain.add_block(gen_complexMap('cSample', concsSample))

    cmain.add_block(gen_calculate(BGElist, sampleList))

    cmain.add_block(gen_check_BGE(expected.pop(0)))

    while expected:
        cmain.add_block(gen_check_eigenzone(expected.pop(0)))

    cmain.add_block(CBlock.make(['return EXIT_SUCCESS;']))

    prog = CProgram()
    prog.add_include('<cstdlib>')
    prog.add_include('\"barsarkagang_tests.h\"')
    prog.add_prelude('using namespace ECHMET;')
    prog.add_prelude('using namespace ECHMET::Barsarkagang;')
    prog.add_function(cmain)

    print(prog)

    try:
        ofh = open(outfile, 'w')
        ofh.write(str(prog))
    except IOError as ex:
        print('Cannot write output: {}'. str(ex))


if __name__ == "__main__":
    parser = make_argparser()
    args = parser.parse_args()

    args_ok = True

    def print_error(s):
        global args_ok
        print('Invalid parameters: {}'.format(s))
        args_ok = False

    if args.input is None:
        print_error('No input file specified')
    if args.output is None:
        print_error('No output file specified')
    if args.generator_path is None:
        print_error('No path to generator executable')
    if args.ECL_path is None:
        print_error('No path to ECHMETCoreLibs binaries')
    if args.LEMNG_path is None:
        print_error('No path to LEMNG binaries')

    if not args_ok:
        sys.exit(1)

    main(args.output, args.input, args.generator_path, args.ECL_path,
         args.LEMNG_path, args.debhue, args.onsfuo, args.viscos)
