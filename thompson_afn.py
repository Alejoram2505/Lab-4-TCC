import sys
import re
from collections import defaultdict, deque

# Utilidades de REGEX -> Postfija

# Operadores y su precedencia
PREC = {
    '*': 3, '+': 3, '?': 3,
    '.': 2,              
    '|': 1
}
RIGHT_ASSOC = {'*', '+', '?'} 

def es_simbolo(c):
    # símbolo "de alfabeto": letra, dígito o epsilon explícita 'ε'
    return (
        c == 'ε' or
        c.isalnum() or
        c in ['_', '#']
    )

def insertar_concatenacion(regex):
    """Inserta el operador '.' donde la concatenación es implícita."""
    res = []
    for i, c1 in enumerate(regex):
        if c1 == ' ':
            continue
        res.append(c1)
        if i == len(regex)-1:
            break
        if regex[i+1] == ' ':
            continue
        c2 = regex[i+1]
        # condiciones donde hay concatenación entre c1 y c2
        if (
            (es_simbolo(c1) or c1 in [')', '*', '+', '?']) and
            (es_simbolo(c2) or c2 in ['(', 'ε'])
        ):
            res.append('.')
    return ''.join(res)

def a_postfija(regex):
    """Shunting-yard."""
    output = []
    opstack = []
    for c in regex:
        if c == ' ':
            continue
        if es_simbolo(c):
            output.append(c)
        elif c == '(':
            opstack.append(c)
        elif c == ')':
            while opstack and opstack[-1] != '(':
                output.append(opstack.pop())
            if not opstack:
                raise ValueError("Paréntesis desbalanceados")
            opstack.pop() 
        elif c in PREC:
            while (opstack and opstack[-1] != '(' and
                   (PREC[opstack[-1]] > PREC[c] or
                    (PREC[opstack[-1]] == PREC[c] and c not in RIGHT_ASSOC))):
                output.append(opstack.pop())
            opstack.append(c)
        else:
            raise ValueError(f"Símbolo no soportado: {c!r}")

    while opstack:
        top = opstack.pop()
        if top in '()':
            raise ValueError("Paréntesis desbalanceados al final")
        output.append(top)
    return ''.join(output)

# Thompson

class Estado:
    __slots__ = ("id", "trans", "eps")
    def __init__(self, id_):
        self.id = id_
        self.trans = defaultdict(set)  
        self.eps = set()               

class AFN:
    def __init__(self, inicio, aceptacion, estados):
        self.inicio = inicio
        self.aceptacion = aceptacion
        self.estados = estados  
def thompson_desde_postfija(post):
    """
    Devuelve AFN usando el algoritmo de Thompson.
    Cada fragmento es (inicio, aceptación)
    """
    next_id = [0]
    def nuevo_estado():
        e = Estado(next_id[0])
        next_id[0] += 1
        return e

    pila = []
    estados = []

    def frag_simbolo(a):
        s = nuevo_estado()
        f = nuevo_estado()
        s.trans[a].add(f)
        estados.extend([s, f])
        pila.append((s, f))

    for c in post:
        if es_simbolo(c):
            # símbolo normal (incluye 'ε' como símbolo de consumo cero)
            if c == 'ε':
                # Thompson para epsilon: s --ε--> f
                s = nuevo_estado()
                f = nuevo_estado()
                s.eps.add(f)
                estados.extend([s, f])
                pila.append((s, f))
            else:
                frag_simbolo(c)
        elif c == '.':
            # concatenación: AB  =>  (a1,a2)·(b1,b2)  une a2 --ε--> b1
            b1, b2 = pila.pop()
            a1, a2 = pila.pop()
            a2.eps.add(b1)
            pila.append((a1, b2))
        elif c == '|':
            # alternancia: A|B, nuevo s,f con ε a ambos inicios y desde ambas acept. a f
            b1, b2 = pila.pop()
            a1, a2 = pila.pop()
            s = nuevo_estado()
            f = nuevo_estado()
            s.eps.update([a1, b1])
            a2.eps.add(f)
            b2.eps.add(f)
            estados.extend([s, f])
            pila.append((s, f))
        elif c == '*':
            # cierre de Kleene
            a1, a2 = pila.pop()
            s = nuevo_estado()
            f = nuevo_estado()
            s.eps.update([a1, f])
            a2.eps.update([a1, f])
            estados.extend([s, f])
            pila.append((s, f))
        elif c == '+':
            # una o más: A+ = A A*
            a1, a2 = pila.pop()
            s = nuevo_estado()
            f = nuevo_estado()
            s.eps.add(a1)
            a2.eps.update([a1, f])
            estados.extend([s, f])
            pila.append((s, f))
        elif c == '?':
            # cero o una: A? = A | ε
            a1, a2 = pila.pop()
            s = nuevo_estado()
            f = nuevo_estado()
            s.eps.update([a1, f])
            a2.eps.add(f)
            estados.extend([s, f])
            pila.append((s, f))
        else:
            raise ValueError(f"Operador inesperado {c}")

    if len(pila) != 1:
        raise ValueError("Expresión inválida (sobran fragmentos)")

    inicio, acept = pila.pop()
    vistos = set()
    q = deque([inicio])
    alcance = []
    while q:
        u = q.popleft()
        if u in vistos:
            continue
        vistos.add(u)
        alcance.append(u)
        for v in u.eps:
            if v not in vistos: q.append(v)
        for _, dests in u.trans.items():
            for v in dests:
                if v not in vistos: q.append(v)

    return AFN(inicio, acept, alcance)

# Simulación de AFN

def epsilon_cierre(estados):
    """Cierre-ε de un conjunto de estados."""
    pila = list(estados)
    cierre = set(estados)
    while pila:
        u = pila.pop()
        for v in u.eps:
            if v not in cierre:
                cierre.add(v)
                pila.append(v)
    return cierre

def mover(estados, simbolo):
    res = set()
    for u in estados:
        for v in u.trans.get(simbolo, ()):
            res.add(v)
    return res

def acepta(afn, cadena):
    actual = epsilon_cierre({afn.inicio})
    for c in cadena:
        actual = epsilon_cierre(mover(actual, c))
        if not actual:
            break
    return afn.aceptacion in actual

# Dibujo del AFN

def dibujar_afn(afn, nombre_png, nombre_dot):

    try:
        import networkx as nx
        import matplotlib.pyplot as plt

        G = nx.DiGraph()
        labels = {}
        for s in afn.estados:
            G.add_node(s.id)
        # recopilar etiquetas por arista (puede haber varias con distinto símbolo)
        edge_labels = defaultdict(list)
        for s in afn.estados:
            for v in s.eps:
                edge_labels[(s.id, v.id)].append('ε')
            for sym, dests in s.trans.items():
                for v in dests:
                    edge_labels[(s.id, v.id)].append(sym)
        for (u,v), labs in edge_labels.items():
            G.add_edge(u, v, label='|'.join(sorted(set(labs))))

        pos = nx.spring_layout(G, seed=42)
        # nodos
        accept_nodes = [afn.aceptacion.id]
        start_node = afn.inicio.id
        others = [n for n in G.nodes() if n not in accept_nodes and n != start_node]

        nx.draw_networkx_nodes(G, pos, nodelist=others)
        nx.draw_networkx_nodes(G, pos, nodelist=[start_node], node_shape='s')
        nx.draw_networkx_nodes(G, pos, nodelist=accept_nodes, node_size=900)
        nx.draw_networkx_labels(G, pos, labels={n:str(n) for n in G.nodes()})
        nx.draw_networkx_edges(G, pos, arrows=True)
        nx.draw_networkx_edge_labels(G, pos,
                                     edge_labels={(u,v):d['label'] for u,v,d in G.edges(data=True)})
        plt.axis('off')
        plt.tight_layout()
        plt.savefig(nombre_png, dpi=170)
        plt.close()
        return f"Imagen generada: {nombre_png}"
    except Exception as e:
        # Fallback: DOT
        try:
            with open(nombre_dot, 'w', encoding='utf-8') as f:
                f.write("digraph NFA {\n  rankdir=LR;\n")
                f.write(f'  node [shape=circle];\n')
                for s in afn.estados:
                    shape = "doublecircle" if s is afn.aceptacion else "circle"
                    f.write(f'  {s.id} [shape={shape}];\n')
                # flecha de inicio
                f.write('  start [shape=point];\n')
                f.write(f'  start -> {afn.inicio.id};\n')
                def w(u, v, lab):
                    return f'  {u} -> {v} [label="{lab}"];\n'
                for s in afn.estados:
                    for v in s.eps:
                        f.write(w(s.id, v.id, 'ε'))
                    for sym, dests in s.trans.items():
                        for v in dests:
                            f.write(w(s.id, v.id, sym))
                f.write("}\n")
            return f"No se pudo dibujar con networkx/matplotlib ({type(e).__name__}). Se generó DOT: {nombre_dot}"
        except Exception as e2:
            return f"No se pudo generar imagen ni DOT: {e2}"

# Entrada / salida

def procesar_linea(linea, idx):
    """Devuelve (regex, w|None)."""
    if ';' in linea:
        izq, der = linea.split(';', 1)
        return izq.strip(), der.strip()
    return linea.strip(), None

def main():
    if len(sys.argv) != 2:
        print("Uso: python thompson_afn.py archivo.txt")
        sys.exit(1)

    ruta = sys.argv[1]
    with open(ruta, 'r', encoding='utf-8') as f:
        lineas = [l.strip() for l in f if l.strip() and not l.strip().startswith('#')]

    for i, linea in enumerate(lineas, start=1):
        regex, w = procesar_linea(linea, i)
        print(f"\n[{i}] r = {regex}")
        if w is None:
            try:
                w = input("  Cadena w (ENTER para vacío): ")
            except EOFError:
                w = ""
        print(f"  w = {w!r}")

        try:
            regl = insertar_concatenacion(regex)
            post = a_postfija(regl)
            afn = thompson_desde_postfija(post)

            # Dibujo
            msg = dibujar_afn(afn, f"nfa_{i}.png", f"nfa_{i}.dot")
            print(" ", msg)

            # Simulación
            ok = acepta(afn, w)
            print("  Resultado:", "sí" if ok else "no")
        except Exception as e:
            print("  Error al procesar la expresión:", e)

if __name__ == "__main__":
    main()
