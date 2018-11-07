#! /usr/bin/env python3

import enum
import json
import os
import subprocess
import sys

import binascii

def _list_to_str(array, formatter, last_formatter=None):
    s = ''

    if len(array) > 0:
        for idx in range(0, len(array) - 1):
            s += formatter(array[idx])
        if last_formatter is None:
            s += formatter(array[-1])
        else:
            s += last_formatter(array[-1])

    return s


def _list_walk(l, handler, handler_last=None):
    if len(l) > 0:
        for idx in range(0, len(l) - 1):
            handler(l[idx])
        if handler_last is None:
            handler(l[-1])
        else:
            handler_last(l[-1])


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
    def __init__(self, lines):
        self._lines = lines

    def add_line(self, line):
        self._lines.append(line)

    @staticmethod
    def _mk_tabs(indent):
        tabs = ''
        for _ in range(0, indent):
            tabs += '\t'

        return tabs

    def __str__(self, indent=0):
        tabs = self._mk_tabs(indent)

        def fmtr(item):
            return '{}{}\n'.format(tabs, item)

        def fmtr_last(item):
            return '{}{}'.format(tabs, item)

        return _list_to_str(self._lines, fmtr, fmtr_last)


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
    s = 'mkRealVec({ { '

    def fmtr(arg):
        return '{}, '.format(arg)

    def last_fmtr(arg):
        return '{}'.format(arg)

    s += _list_to_str(items, fmtr, last_fmtr)

    s += ' } })'

    return s


def gen_complexes(cpxs):
    if cpxs is None:
        return 'nullptr'

    if len(cpxs) == 0:
        return 'noComplexes()'


def gen_constituent(name, ctype, chargeLow, chargeHigh, pKas, mobilities, viscosity,
                    complexations):
    ctuent = StructInitializer('SysComp::InConstituent', cname_from_name(name))

    ctuent.add_item(gen_ctype(ctype))
    ctuent.add_item(gen_fixedString(name))
    ctuent.add_item(chargeLow)
    ctuent.add_item(chargeHigh)
    ctuent.add_item(gen_realvec(pKas))
    ctuent.add_item(gen_realvec(mobilities))
    ctuent.add_item(gen_complexes(complexations))
    ctuent.add_item(viscosity)

    return ctuent


def gen_complexMap(name, m):
    block = CBlock(['CMapping {} = {{'.format(name)])

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
    block = CBlock(['const auto r = calculate('])

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
    return CBlock(['checkBGE(r, {}, {}, {}, {});'.format(
        expected[0], expected[1], expected[2], expected[3])])


def gen_check_eigenzone(expected):
    return CBlock(['checkEigenzone(r.eigenzone, {}, {}, {}, {});'.format(
        expected[0], expected[1], expected[2], expected[3])])


def process_complexes(root):
    return []


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
                return process_complexes(c['complexForms'])
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
    if is_corr > 0:
        params.extend(['1', '1', '1'])
    else:
        params.extend(['0', '0', '0'])

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


def main(infile, genpath, ecl_path, lemng_path, is_corr):
    is_corr = int(is_corr)

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

    while len(expected) > 0:
        cmain.add_block(gen_check_eigenzone(expected.pop(0)))

    cmain.add_block(CBlock(['return EXIT_SUCCESS;']))

    prog = CProgram()
    prog.add_include('<cstdlib>')
    prog.add_include('\"barsarkagang_tests.h\"')
    prog.add_prelude('using namespace ECHMET;')
    prog.add_prelude('using namespace ECHMET::Barsarkagang;')
    prog.add_function(cmain)

    print(prog)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
