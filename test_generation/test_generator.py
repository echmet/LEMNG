#! /usr/bin/env python3

import argparse
import enum
import json
import os
import subprocess
import sys
from enum import Enum


INDENT_SPACER = '\t'


class RefDataResult(Enum):
    Normal = 1,
    Fail = 2,
    Oscillating = 3


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

        last = len(self._code_blocks) - 1
        for i, block in enumerate(self._code_blocks):
            s += '{}\n{}'.format(block, '' if i == last else '\n')

        s += '}'

        return s

    def add_block(self, block):
        self._code_blocks.append(block.__str__(1))


class CBlock:
    @staticmethod
    def make(items):
        block = CBlock()
        for i in items:
            block.add_item(i)
        return block

    def __init__(self):
        self._elems = []

    def add_item(self, elem):
        self._elems.append(elem)

    @staticmethod
    def _mk_spcs(indent):
        spcs = ''
        for _ in range(0, indent):
            spcs += INDENT_SPACER

        return spcs

    def __str__(self, indent=0):
        spcs = self._mk_spcs(indent)

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
        spcs = self._mk_spcs(indent)

        s = '{}{} {}{{\n'.format(spcs, self._stype, self._name)

        def fmtr(code):
            return '{}{}{},\n'.format(INDENT_SPACER, spcs, code)

        def last_fmtr(code):
            return '{}{}{}\n'.format(INDENT_SPACER, spcs, code)

        s += _list_to_str(self._items, fmtr, last_fmtr)

        s += '{}}};'.format(spcs)

        return s

    def add_item(self, code):
        self._items.append(code)


def cname_from_name(name):
    return name.lower().replace(' ', '_').replace('-', '__')


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
    s = 'mkRealVec( '
    s += _list_to_initer_list(items)
    s += ' )'

    return s


def gen_complexes(cpxs, name):
    if not cpxs:
        return None

    genfunc_name = 'gen_complexforms_' + cname_from_name(name)
    genfunc = CFunction('SysComp::InCFVec *', genfunc_name, [])
    outerblock = CBlock().make(['const ComplexDef cDef = {'])

    incf_last = len(cpxs) - 1
    for i, incf in enumerate(cpxs):
        ocfblock = CBlock.make(['{ /* InComplexForm c-tor begin */'])

        cfblock = CBlock()
        cfblock.add_item('{},'.format(incf['nucleusCharge']))

        cfblock.add_item('/* InLGVec */')
        cfblock.add_item('{')

        inlg_last = len(incf['ligandGroups']) - 1
        for j, inlg in enumerate(incf['ligandGroups']):
            incfblock = CBlock()

            incfblock.add_item('{ /* InLigandGroup c-tor begin */')

            lgblock = CBlock()

            lgblock.add_item('/* InLFVec */')
            lgblock.add_item('{')

            inlf_last = len(inlg['ligands']) - 1
            for k, inlf in enumerate(inlg['ligands']):
                lfblock = CBlock()

                lfblock.add_item('{ /* InLigandForm c-tor begin */')

                inlfblock = CBlock()

                inlfblock.add_item('\"{}\",'.format(inlf['name']))
                inlfblock.add_item('{},'.format(inlf['charge']))
                inlfblock.add_item('{},'.format(inlf['maxCount']))
                inlfblock.add_item('{},'.format(_list_to_initer_list(inlf['pBs'])))
                inlfblock.add_item('{}'.format(_list_to_initer_list(inlf['mobilities'])))

                lfblock.add_item(inlfblock)

                if k != inlf_last:
                    lfblock.add_item('}, /* InLigandForm c-tor end */')
                else:
                    lfblock.add_item('} /* InLigandForm c-tor end */')

                lgblock.add_item(lfblock)

            lgblock.add_item('}')
            incfblock.add_item(lgblock)

            if j != inlg_last:
                incfblock.add_item('}, /* InLigandGroup c-tor end */')
            else:
                incfblock.add_item('} /* InLigandGroup c-tor end */')

            cfblock.add_item(incfblock)

        cfblock.add_item('}')

        ocfblock.add_item(cfblock)
        ocfblock.add_item('}}{} /* InComplexForm c-tor end */'.format('' if i == incf_last else ','))

        outerblock.add_item(ocfblock)

    outerblock.add_item('};')

    genfunc.add_block(outerblock)
    genfunc.add_block(CBlock.make(['return buildComplexes(cDef);']))

    return (genfunc_name, genfunc)


def gen_constituent(name, ctype, chargeLow, chargeHigh, pKas, mobilities,
                    viscosity, complexations):
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


def gen_concMap(name, m):
    block = CBlock.make(['CMapping {} = {{'.format(name)])

    e = enumerate(m)

    for idx, k in e:
        v = m[k]
        if idx == len(m) - 1:
            block.add_item('{}{{ \"{}\", {} }}'.format(INDENT_SPACER, k, v))
        else:
            block.add_item('{}{{ \"{}\", {} }},'.format(INDENT_SPACER, k, v))

    block.add_item('};')

    return block


def gen_calculate(BGElist, sampleList, is_corr, should_oscillate):
    block = CBlock.make(['const auto r = calculate('])

    def mk_vec(l):
        block.add_item(INDENT_SPACER + '{')

        def hndl(item):
            block.add_item(INDENT_SPACER + INDENT_SPACER + '{},'.format(cname_from_name(item)))

        def hndl_last(item):
            block.add_item(INDENT_SPACER + INDENT_SPACER + '{}'.format(cname_from_name(item)))

        _list_walk(l, hndl, hndl_last)

        block.add_item(INDENT_SPACER + '},')

    mk_vec(BGElist)
    mk_vec(sampleList)

    block.add_item(INDENT_SPACER + 'cBGE, cSample,')

    def corr_set(n):
        return 'true' if is_corr & n else 'false'

    block.add_item(INDENT_SPACER + '{0}, {1}, {2}, {3});'.format(
                   corr_set(1), corr_set(2), corr_set(4),
                   'true' if should_oscillate else 'false'))

    return block


def gen_check_BGE(expected):
    return CBlock.make(['checkBGE(r, {}, {}, {}, {});'.format(
        expected[0], expected[1], expected[2], expected[3])])


def gen_check_eigenzone(ez_id, expected):
    return CBlock.make(['checkEigenzone({}, r.eigenzones, {}, {}, {}, {}, {});'.format(
        ez_id, expected[0], expected[1], expected[2], expected[3], expected[4])])


def process_input(root):
    ctuents = root['constituents']

    c_ctuents = []
    concsBGE = dict()
    concsSample = dict()
    BGElist = []
    sampleList = []
    complex_generators = dict()

    for c in ctuents:
        name = c['name']

        if 'complexForms' in c:
            complex_generators[name] = gen_complexes(c['complexForms'], name)

        def mk_cpx():
            if name in complex_generators:
                cg = complex_generators[name]
                if cg is None:
                    return 'noComplexes()'
                return complex_generators[name][0] + '()'
            return 'nullptr'

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

        concsSample[c['name']] = c['concentrationSample']
        sampleList.append(c['name'])
        if c['role'] == 'B':
            BGElist.append(c['name'])
            concsBGE[c['name']] = c['concentrationBGE']

    return (c_ctuents, concsBGE, concsSample, BGElist, sampleList,
            complex_generators)


def results_file(infile):
    los = infile.rfind('/')
    lod = infile.rfind('.')

    return infile[los+1:lod] + '_results.txt'


def get_expected_results(genpath, infile, ecl_path, lemng_path, is_corr, silent):
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
            1) mobility 2) uEMD 3) a2t 4) pH 5) Conductivity
        """

        data = []
        while len(lines) > 0:
            ez = []
            for _ in range(0, 5):
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

    fhout = open(os.devnull, 'w') if silent else None

    try:
        ret = subprocess.run(params, stdout=fhout, stderr=fhout)
        if ret.returncode == 0x66:
            return (RefDataResult.Oscillating, [])
        elif ret.returncode != 0:
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

    return (RefDataResult.Normal, data)


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
    parser.add_argument('--silent', help='Do not display any output from reference generator', action='store_true')

    parser.set_defaults(debhue=False)
    parser.set_defaults(onsfuo=False)
    parser.set_defaults(viscos=False)
    parser.set_defaults(silent=False)

    return parser


def print_cmake_line(outfile):
    spcs = '    '

    fidx = outfile.rfind('/') + 1
    tidx = outfile.rfind('.')

    if tidx <= fidx:
        tidx = len(outfile)

    tag = outfile[fidx:tidx]
    basename = outfile[fidx:]

    print(spcs + 'add_executable({}_exe src/tests/{})'.format(tag, basename))

    indent = ''
    for _ in range(0, len('target_link_libraries({}_exe '.format(tag))):
        indent += ' '

    print(spcs + 'target_link_libraries({}_exe PRIVATE LEMNG'.format(tag))
    print(spcs + indent + 'PRIVATE ECHMETShared')
    print(spcs + indent + 'PRIVATE SysComp)')

    print(spcs + 'add_test({0} {0}_exe)\n'.format(tag))


def main(outfile, infile, genpath, ecl_path, lemng_path, debhue, onsfuo,
         viscos, silent):
    is_corr = 0
    if debhue:
        is_corr += 1
    if onsfuo:
        is_corr += 2
    if viscos:
        is_corr += 4

    refres, expected = get_expected_results(genpath, infile, ecl_path,
                                            lemng_path, is_corr, silent)

    fh = open(infile, 'r')
    sysdef = json.load(fh)

    (c_list, concsBGE, concsSample, BGElist,
     sampleList, complex_generators) = process_input(sysdef)

    prog = CProgram()
    prog.add_include('<cstdlib>')
    prog.add_include('\"barsarkagang_tests.h\"')
    prog.add_prelude('using namespace ECHMET;')
    prog.add_prelude('using namespace ECHMET::Barsarkagang;')

    cmain = CFunction('int', 'main', [CFunctionArg('int', ''),
                                      CFunctionArg('char **', '')])

    for c in c_list:
        cmain.add_block(c)

    cmain.add_block(gen_concMap('cBGE', concsBGE))
    cmain.add_block(gen_concMap('cSample', concsSample))

    if refres == RefDataResult.Oscillating:
        cmain.add_block(gen_calculate(BGElist, sampleList, is_corr, True))
        cmain.add_block(CBlock.make(['(void)r;']))
    else:
        cmain.add_block(gen_calculate(BGElist, sampleList, is_corr, False))

        cmain.add_block(gen_check_BGE(expected.pop(0)))

        ez_id = 1
        while expected:
            cmain.add_block(gen_check_eigenzone(ez_id, expected.pop(0)))
            ez_id += 1

    cmain.add_block(CBlock.make(['return EXIT_SUCCESS;']))

    for _, gf in complex_generators.items():
        if gf is not None:
            prog.add_function(gf[1])

    prog.add_function(cmain)

    try:
        ofh = open(outfile, 'w')
        ofh.write(str(prog))
    except IOError as ex:
        print('Cannot write output: {}'.format(ex))
        return

    print_cmake_line(outfile)


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
         args.LEMNG_path, args.debhue, args.onsfuo, args.viscos, args.silent)
