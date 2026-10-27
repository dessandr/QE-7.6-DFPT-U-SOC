"""Resolve the qe-pert SOC+U preprocessor macros in a Fortran source file.

The q-operator path (__SOC_U_DFPT_Q_OPERATOR_VALIDATE) is the ported feature,
so it is treated as defined; the control/diagnostic twins and staged guard
modes are treated as undefined. Directives that test any other macro
(__MPI, __CUDA, __PH_ACC, ...) are kept verbatim, while known macros nested
inside them are still resolved. Usage:
    python3 -I soc_u_port_resolve_macros.py IN_FILE OUT_FILE
"""
import re
import sys

DEFINED = {'__SOC_U_DFPT_Q_OPERATOR_VALIDATE'}
UNDEFINED = {
    '__SOC_U_DFPT_FIXED_CG_VALIDATE', '__SOC_U_DFPT_LEGACY_MIX_VALIDATE',
    '__SOC_U_DFPT_ORIGINAL_CONV_VALIDATE', '__SOC_U_DFPT_JOINT_HXC_ONLY_VALIDATE',
    '__SOC_U_DFPT_NOHUB_FEEDBACK_VALIDATE', '__SOC_U_PW_RESIDUAL_DEBUG',
    '__SOC_U_DFPT_Q_BARE_VALIDATE', '__SOC_U_DFPT_FULL_VALIDATE',
    '__SOC_U_DFPT_BARE_VALIDATE',
}
KNOWN = DEFINED | UNDEFINED
DIRECTIVE = re.compile(r'^\s*#\s*(if|ifdef|ifndef|elif|else|endif)\b(.*)$')
MACRO = re.compile(r'\b(__\w+)\b')


def evaluate(kind, expr):
    """True/False if the condition uses only known macros, else None."""
    if kind in ('ifdef', 'ifndef'):
        name = expr.strip()
        if name not in KNOWN:
            return None
        value = name in DEFINED
        return value if kind == 'ifdef' else not value
    names = set(MACRO.findall(expr))
    if not names or not names <= KNOWN:
        return None
    py = re.sub(r'defined\s*\(\s*(\w+)\s*\)', lambda m: str(m.group(1) in DEFINED), expr)
    py = re.sub(r'defined\s+(\w+)', lambda m: str(m.group(1) in DEFINED), py)
    py = py.replace('&&', ' and ').replace('||', ' or ').replace('!', ' not ')
    if re.search(r'[^\sA-Za-z()]', py) or set(re.findall(r'[A-Za-z]+', py)) - {'True', 'False', 'and', 'or', 'not'}:
        return None
    return bool(eval(py, {'__builtins__': {}}, {}))


def resolve(lines):
    out = []
    # Each frame: [resolved, emitting_parent, branch_taken, current_branch_on]
    stack = []

    def active():
        return all(f[3] for f in stack)

    for line in lines:
        m = DIRECTIVE.match(line)
        if not m:
            if active():
                out.append(line)
            continue
        kind, rest = m.group(1), m.group(2)
        if kind in ('if', 'ifdef', 'ifndef'):
            value = evaluate(kind, rest)
            if value is None:
                stack.append([False, active(), None, True])
                if active():
                    out.append(line)
            else:
                stack.append([True, active(), value, value])
        elif kind == 'elif':
            frame = stack[-1]
            value = evaluate('if', rest)
            if not frame[0]:
                if value is not None:
                    raise ValueError(f'known macro in #elif of unknown #if: {line!r}')
                if all(f[3] for f in stack[:-1]):
                    out.append(line)
            else:
                if value is None:
                    raise ValueError(f'unknown macro in #elif of known #if: {line!r}')
                frame[3] = (not frame[2]) and value
                frame[2] = frame[2] or value
        elif kind == 'else':
            frame = stack[-1]
            if not frame[0]:
                if all(f[3] for f in stack[:-1]):
                    out.append(line)
            else:
                frame[3] = not frame[2]
                frame[2] = True
        elif kind == 'endif':
            frame = stack.pop()
            if not frame[0] and active():
                out.append(line)
    if stack:
        raise ValueError('unbalanced preprocessor directives')
    return out


def main():
    src, dst = sys.argv[1:3]
    with open(src) as f:
        lines = f.readlines()
    with open(dst, 'w') as f:
        f.writelines(resolve(lines))


if __name__ == '__main__':
    main()
