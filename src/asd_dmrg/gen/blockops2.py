from __future__ import print_function
import sys
import math as m

class Op:
    spin = False
    creation = False

    def __init__(self, sp, cr):
        self.spin = sp
        self.creation = cr

    def conjugate(self):
        return Op(self.spin, not self.creation)

    def GammaSQ(self):
        out = "GammaSQ::"
        out += "Create" if self.creation else "Annihilate"
        out += "Alpha" if self.spin else "Beta"
        return out

    def __str__(self):
        label = "A" if self.spin else "B"
        if (self.creation):
          label += "^t"
        return label

    def __eq__(self, other):
        return self.spin==other.spin and self.creation==other.creation

    def __ne__(self, other):
        return not self==other

class OpList:
    ops = []

    def __init__(self, oplist):
        self.ops = oplist

    def conjugate(self):
        return OpList([k.conjugate() for k in reversed(self.ops)])

    def reverse(self):
        return OpList([k for k in reversed(self.ops)])

    def __str__(self):
        out = ""
        if (len(self.ops)!=0):
            for o in self.ops:
                out += " %s" % o
            out += " "
            return out
        else:
            return "  I"

    def __len__(self):
        return len(self.ops)

    def __eq__(self, other):
        if len(self)!=len(other):
            return False
        else:
            for s, o in zip(self.ops, other.ops):
                if (s!=o):
                    return False
            return True

    def __getitem__(self, i):
        return self.ops[i]

    def split(self, partition):
        leftops = []
        rightops = []
        for i, o in enumerate(self.ops):
              if ((partition & (1 << i))):
                  leftops.append(o)
              else:
                  rightops.append(o)
        return (OpList(leftops), OpList(rightops))

    def alpha_change(self):
        out = 0
        for o in self.ops:
            if o.spin:
                if o.creation:
                    out += 1
                else:
                    out -= 1
        if out==0:
            return "  "
        elif out > 0:
            return "+%s" % out
        else:
            return "%2s" % out

    def beta_change(self):
        out = 0
        for o in self.ops:
            if not o.spin:
                if o.creation:
                    out += 1
                else:
                    out -= 1
        if out==0:
            return "  "
        elif out > 0:
            return "+%s" % out
        else:
            return "%2s" % out

def CreaA():
    return Op(True, True)

def AnniA():
    return Op(True, False)

def CreaB():
    return Op(False, True)

def AnniB():
    return Op(False, False)

allowed_list = [
    OpList([CreaA()]),
    OpList([CreaB()]),
    OpList([AnniA(), AnniA()]),
    OpList([AnniB(), AnniB()]),
    OpList([AnniB(), AnniA()]),
    OpList([CreaA(), AnniA()]),
    OpList([CreaB(), AnniB()]),
    OpList([CreaB(), AnniA()]),
    OpList([CreaB(), CreaA(), AnniA()]),
    OpList([CreaA(), CreaA(), AnniA()]),
    OpList([CreaB(), CreaB(), AnniB()]),
    OpList([CreaA(), CreaB(), AnniB()])
]

class Term(OpList):
    rev = False
    conj = False

    def __init__(self, oplist):
        self.ops = oplist
        if len(oplist)==0:
            self.rev, self.conj = False, False
        elif oplist in allowed_list:
            self.rev, self.conj = False, False
        elif oplist.conjugate() in allowed_list:
            self.rev, self.conj = False, True
        elif oplist.reverse() in allowed_list:
            self.rev, self.conj = True, False
        elif oplist.reverse().conjugate() in allowed_list:
            self.rev, self.conj = True, True
        else:
            raise Exception("This should not have happened!")

    def get_allowed(self):
        if not self.rev and not self.conj:
            return self.ops
        elif self.rev and not self.conj:
            return self.reverse()
        elif not self.rev and self.conj:
            return self.conjugate()
        else:
            return self.reverse().conjugate()

    def factor(self):
        return 1.0 if not self.rev else -1.0

    def matel(self, label):
        bra = label + "'"
        ket = label
        out = "<%s|%s|%s>" % (bra, OpList.__str__(self), ket)
        return out

    def GammaSQ(self):
        oplist = self.get_allowed()
        out = "%s" % oplist[0].GammaSQ()
        for i in range(1, len(oplist)):
            out += ", %s" % oplist[i].GammaSQ()
        return out

    def __str__(self):
        out = "-" if self.rev else ""
        out += "%s" % OpList.__str__(self)
        return str(out)

    def __len__(self):
        return len(self.ops)

level = 0
def indent():
    return "  "*level

def open_code_block():
    global level
    level += 1

def close_code_block():
    global level
    level -= 1
    print("%s}" % indent())

input_indices = [ "i", "j", "k" ]
left_indices = [ "a", "b", "c" ]
right_indices = [ "p", "q", "r" ]

norb = "norb"
lnorb = "lnorb"
rnorb = "rnorb"
loffset = "loffset"
roffset = "roffset"
mo2e = "mo2e"

class ContractedOperator(OpList):
    integrals = ""
    factor = 1.0

    def __init__(self, ops, integrals, factor):
        self.ops = ops
        self.integrals = integrals
        self.factor = factor

class Integrals:
    partition = bin(0)
    size = 1
    reverseL = False
    reverseR = False
    int_type = "J"

    def index(self, indices):
        out = ""
        for i, index in enumerate(indices):
            out += "%s%s" % ("" if i==0 else ", ", index)
        return out

    def __init__(self, partition, size, rL, rR, it):
        self.size = size
        self.partition = bin(partition)[2:]
        self.partition = "0"*(size - len(self.partition)) + self.partition
        self.reverseL = rL
        self.reverseR = rR
        self.int_type = it

    def __call__(self, ptr_name, Cindices, Lindices, Rindices):
        if self.reverseL:
            Lindices.reverse()
        if self.reverseR:
            Rindices.reverse()

        size = len(Lindices) + len(Rindices)

        if len(self.partition) != size:
            self.partition = "0"*(size - len(self.partition)) + self.partition

        total_indices = []
        for i in reversed(range(size)):
            if self.partition[i] == "1":
                total_indices.append(Lindices.pop(0) + "+" + loffset)
            else:
                total_indices.append(Rindices.pop(0) + "+" + roffset)

        if self.int_type == "QJ":
            assert size == 2
            out = "%s(%s)" % (ptr_name, self.index( [total_indices[0], Cindices[0], total_indices[1], Cindices[1]] ))
        elif self.int_type == "QK":
            assert size == 2
            out = "%s(%s)" % (ptr_name, self.index( [total_indices[0], total_indices[1], Cindices[1], Cindices[0]] ))
        elif self.int_type == "QJK":
            assert size == 2
            out = "(%s(%s) - %s(%s))" % (ptr_name, self.index([total_indices[0], Cindices[0], total_indices[1], Cindices[1]]),
                                         ptr_name, self.index([Cindices[0], total_indices[0], total_indices[1], Cindices[1]]))
        elif self.int_type == "P":
            assert size == 2
            out = "%s(%s)" % (ptr_name, self.index([total_indices[0], total_indices[1], Cindices[1], Cindices[0]]))
        elif self.int_type == "PJK":
            assert size == 2
            out = "(%s(%s) - %s(%s))" % (ptr_name, self.index([total_indices[0], total_indices[1], Cindices[1], Cindices[0]]),
                                         ptr_name, self.index([total_indices[0], total_indices[1], Cindices[0], Cindices[1]]))
        elif self.int_type == "S":
            assert size == 3
            out = "%s(%s)" % (ptr_name, self.index([Cindices[0], total_indices[1], total_indices[0], total_indices[2]]))
        elif self.int_type == "SJK":
            assert size == 3
            out = "(%s(%s) - %s(%s))" % (ptr_name, self.index([Cindices[0], total_indices[1], total_indices[0], total_indices[2]]),
                                         ptr_name, self.index([Cindices[0], total_indices[0], total_indices[1], total_indices[2]]))
        else:
            raise Exception("Undefined integral type!")

        return out


class OperatorProduct:
    left = Term([])
    right = Term([])
    fac = 0.0
    integrals = Integrals(0, 1, False, False, "")

    def __init__(self, lterm, rterm, integrals, fac):
        self.left = lterm
        self.right = rterm
        self.integrals = integrals
        self.fac = fac

    def factor(self):
        return self.fac * self.left.factor() * self.right.factor()

    def action(self):
        return self.left.alpha_change() + self.left.beta_change() + self.right.alpha_change() + self.right.beta_change()

    def compatible(self, other):
        return OpList.__eq__(self.left, other.left) and OpList.__eq__(self.right, other.right) and \
               (self.integrals.int_type == "S" and other.integrals.int_type == "S" or \
                self.integrals.int_type == "P" and other.integrals.int_type == "P")

    def combine(self, other):
        fac = 0.0
        partition = 0
        if self.integrals.int_type == "S":
            if self.integrals.partition in [ "101", "001" ]:
                assert other.integrals.partition in [ "110", "010" ]
                fac = self.factor()
                for i, p in enumerate(reversed(self.integrals.partition)):
                    if p == "1":
                        partition += 1 << i
            elif self.integrals.partition in [ "110", "010" ]:
                assert other.integrals.partition in [ "101", "001" ]
                fac = other.factor()
                for i, p in enumerate(reversed(other.integrals.partition)):
                    if p == "1":
                        partition += 1 << i
            else:
                raise Exception("Incompatible!")
            return OperatorProduct(self.left, self.right, Integrals(partition, 3, False, False, "SJK"), fac)
        elif self.integrals.int_type == "P":
            if self.integrals.partition == "01":
                assert other.integrals.partition == "10"
                fac = self.factor()
                for i, p in enumerate(reversed(self.integrals.partition)):
                    if p == "1":
                        partition += 1 << i
            elif self.integrals.partition == "10":
                assert other.integrals.partition == "01"
                fac = other.factor()
                for i, p in enumerate(reversed(other.integrals.partition)):
                    if p == "1":
                        partition += 1 << i
            else:
                raise Exception("Incompatible!")
            return OperatorProduct(self.left, self.right, Integrals(partition, 2, self.left.conj, self.left.rev, "PJK"), fac)

    def __str__(self):
        return "%1.1f <L'|%s|L> (x) <R'|%s|R>" % (self.factor(), self.left, self.right)

def generate_operator(opname, contracted_operators, ninput):
    def_inp_string = "const int %s" % input_indices[0]
    inp_string = input_indices[0]
    for i in range(1,ninput):
        inp_string += ", %s" % input_indices[i]
        def_inp_string += ", const int %s" % input_indices[i]
    print("%sshared_ptr<BlockSparseMatrix> BlockOperators2::%s(BlockKey bk, %s) const {" % (indent(), opname, def_inp_string))
    open_code_block()

    # action needs to be the same for all of them, so figure out action from the first one
    action = (contracted_operators[0].alpha_change(), contracted_operators[0].beta_change())
    global_action = action[0] + action[1]

    opcollection = {}
    opcollection["    %s" % global_action] = []
    opcollection["%s    " % global_action] = []
    # preprocess into OperatorProducts
    for op in contracted_operators:
        size = len(op)
        npart = 1 << size
        for partitioning in range(1, npart - 1):
            # used to figure out the sign of rearranging the operators so all of the "right" operators are on the right
            split_factor = -1.0 if partitioning==2 or partitioning==5 else 1.0
            (left, right) = op.split(partitioning)
            (lterm, rterm) = (Term(left), Term(right))
            op_prod = OperatorProduct(lterm, rterm, Integrals(partitioning, 4 - ninput, lterm.rev ^ lterm.conj, rterm.rev ^ rterm.conj, op.integrals), op.factor*split_factor)
            if op_prod.action() in opcollection:
                opcollection[op_prod.action()].append(op_prod)
            else:
                opcollection[op_prod.action()] = [ op_prod ]
    diag = action[0]=="  " and action[1]=="  "

    if not diag:
        print("%sBlockKey target_bk(bk.nelea%s,bk.neleb%s);" % (indent(), action[0].strip(), action[1].strip()))
        print("%sassert(blocks_->contains(target_bk));" % indent())
        print()

    svec = "pvec" if diag else "source_pvec"
    tvec = "pvec" if diag else "target_pvec"

    print("%sconst vector<DMRG::BlockPair>& %s = blocks_->blockpairs(bk);" % (indent(), svec))
    if not diag:
        print("%sconst vector<DMRG::BlockPair>& %s = blocks_->blockpairs(target_bk);" % (indent(), tvec))

    print()

    sinfo = "binfo" if diag else "source_info"
    tinfo = "binfo" if diag else "target_info"

    print("%sconst BlockInfo %s = blocks_->blockinfo(bk);" % (indent(), sinfo))
    if not diag:
        print("%sconst BlockInfo %s = blocks_->blockinfo(target_bk);" % (indent(), tinfo))

    print()

    print("%sconst int %s = jop_->nocc();" % (indent(), norb))
    print("%sconst int %s = blocks_->left_block()->norb();" %(indent(), lnorb))
    print("%sconst int %s = blocks_->right_block()->norb();" %(indent(), rnorb))
    print("%sconst int %s = %s - (%s + %s); // convenience variable for offset of left orbitals from zero" % (indent(), loffset, norb, lnorb, rnorb))
    print("%sconst int %s = %s - %s; // convenience variable for offset of right orbitals from zero" % (indent(), roffset, norb, rnorb))
    print()

    print("%sconst btas::TensorView4<double> mo2e = btas::make_view(btas::CRange<4>(%s,%s,%s,%s), jop_->mo2e()->storage());" % (indent(), norb, norb, norb, norb))
    print()

    print("%smap<pair<size_t, size_t>, shared_ptr<Matrix>> out;" % (indent()))
    print()
    print("%sfor (auto& spair : %s) {" % (indent(), svec))
    open_code_block()

    print("%s// phase accumulated by moving an operator past the whole left ket block" % indent())
    print("%sconst int left_phase = 1 - (((spair.left.nelea+spair.left.neleb)%%2) << 1);" % indent());

    for key, collection in opcollection.iteritems():
        pure_left = False
        pure_right = False
        terms_computed = ""
        if key == ("    %s"%global_action):
            pure_right = True
            terms_computed += " I (x) %s" % opname
        if key == ("%s    "%global_action):
            pure_left = True
            if terms_computed != "":
                terms_computed += " + "
            terms_computed += " %s (x) I " % opname
        for i, t in enumerate(collection):
            sgn = ""
            if t.factor() > 0:
                sgn = " + "
            elif t.factor() < 0:
                sgn = " - "
            terms_computed += "%s%1.1f %s (x) %s" % (sgn, abs(t.factor()), t.left.matel("L"), t.right.matel("R"))
        print()
        print("%s{ // %s" % (indent(), terms_computed.strip()))
        open_code_block()

        la = key[0:2]
        lb = key[2:4]
        ra = key[4:6]
        rb = key[6:8]

        #print(la, lb, ra, rb)

        both_diag = (la, lb, ra, rb) == ("  ", "  ", "  ", "  ")

        if not both_diag:
            if (la, lb)==("  ", "  "):
                print("%sconst BlockKey left_target = spair.left.key();" % indent())
            else:
                print("%sconst BlockKey left_target(spair.left.nelea%s, spair.left.neleb%s);" % (indent(), la.strip(), lb.strip()))
            if (ra, rb)==("  ", "  "):
                print("%sconst BlockKey right_target = spair.right.key();" % indent())
            else:
                print("%sconst BlockKey right_target(spair.right.nelea%s, spair.right.neleb%s);" % (indent(), ra.strip(), rb.strip()))
            print()
            print("%sauto iter = find_if(%s.begin(), %s.end(), [&left_target, &right_target] (DMRG::BlockPair bp)" % (indent(), tvec, tvec))
            print("%s  { return make_pair(bp.left.key(), bp.right.key())==make_pair(left_target, right_target); }" % indent())
            print("%s);" % indent())

            print("%sif(iter!=%s.end()) {" % (indent(), tvec))
            open_code_block()

            print("%sDMRG::BlockPair tpair = *iter;" % indent())
            print()

        tpair = "spair" if both_diag else "tpair"
        nops = 0

        print("%sauto out_block = make_shared<Matrix>(%s.nstates(), spair.nstates(), true);" % (indent(), tpair))

        if (pure_right):
            nops += 1
            print()
            print("%s// I (x) %s" % (indent(), opname))
            print("%sMatrix Lident(spair.left.nstates, spair.left.nstates, true); Lident.unit();" % indent())
            print("%sMatrix Rterms = *right_ops_->%s_as_matrix(spair.right.key(), %s);" % (indent(), opname, inp_string))
            print()
            phase = "left_phase" if (ninput%2==1) else "1.0"
            print("%skronecker_product(%s, false, Rterms, false, Lident, *out_block);" % (indent(), phase))

        if (pure_left):
            nops += 1
            print()
            print("%s// %s (x) I" % (indent(), opname))
            print("%sMatrix Lterms = *left_ops_->%s_as_matrix(spair.left.key(), %s);" % (indent(), opname, inp_string))
            print("%sMatrix Rident(spair.right.nstates, spair.right.nstates, true); Rident.unit();" % indent())
            print()
            print("%skronecker_product(1.0, false, Rident, false, Lterms, *out_block);" % (indent()))

        # preprocess collection into terms that can be combined
        combined_collection = []
        while len(collection) != 0:
            seed = collection.pop(0)
            for i in reversed(range(len(collection))):
                if seed.compatible(collection[i]):
                    seed = seed.combine(collection.pop(i))
            combined_collection.append(seed)

        processed_collection = []
        while len(combined_collection) != 0:
            seed =  [ combined_collection.pop(0) ]
            for i in reversed(range(len(combined_collection))):
                #if len(seed[0].left) ==1 and len(combined_collection[i].left) ==1 and len(seed[0].right)==2 and OpList.__eq__(seed[0].left, combined_collection[i].left) or \
                #   len(seed[0].right)==1 and len(combined_collection[i].right)==1 and len(seed[0].left) ==2 and OpList.__eq__(seed[0].right, combined_collection[i].right):
                if len(seed[0].left) ==1 and len(combined_collection[i].left) ==1 and OpList.__eq__(seed[0].left, combined_collection[i].left) or \
                   len(seed[0].right)==1 and len(combined_collection[i].right)==1 and OpList.__eq__(seed[0].right, combined_collection[i].right):
                    seed.append(combined_collection.pop(i))
            processed_collection.append(seed)

        nops += len(processed_collection)

        for oprodvec in processed_collection:
            oprod = oprodvec[0]

            if nops > 1:
                term_string = [ "%s"%x for x in oprodvec ]
                print()
                print("%s{ // %s" % (indent(), term_string))
                open_code_block()

            shorter = "left" if len(oprod.left)==1 and len(oprod.right)==2 else "right"

            lgammas = [ ]
            rgammas = [ ]
            lbra, lket = "tpair", "spair"
            if oprod.left.conj:
                lbra, lket = lket, lbra
            rbra, rket = "tpair", "spair"
            if oprod.right.conj:
                rbra, rket = rket, rbra
            if len(oprodvec)==1:
                lgammas.append( ("Lgamma", oprod.left.GammaSQ(), lbra, lket) )
                rgammas.append( ("Rgamma", oprod.right.GammaSQ(), rbra, rket) )
            elif len(oprodvec)==2:
                if shorter == "left":
                    lgammas.append( ("Lgamma", oprod.left.GammaSQ(), lbra, lket) )
                    rgammas.append( ("Rgamma1", oprodvec[0].right.GammaSQ(), rbra, rket) )
                    rgammas.append( ("Rgamma2", oprodvec[1].right.GammaSQ(), rbra, rket) )
                else:
                    rgammas.append( ("Rgamma", oprod.right.GammaSQ(), rbra, rket) )
                    lgammas.append( ("Lgamma1", oprodvec[0].left.GammaSQ(), lbra, lket) )
                    lgammas.append( ("Lgamma2", oprodvec[1].left.GammaSQ(), lbra, lket) )
            else:
                raise Exception("NO")

            # find individual tensors
            if nops <= 1:
                print()

            for LG in lgammas:
                print("%sshared_ptr<const btas::Tensor3<double>> %s = blocks_->left_block()->coupling({%s}).at({%s.left.key(),%s.left.key()}).data;" % (indent(), LG[0], LG[1], LG[2], LG[3]))
            for RG in rgammas:
                print("%sshared_ptr<const btas::Tensor3<double>> %s = blocks_->right_block()->coupling({%s}).at({%s.right.key(),%s.right.key()}).data;" % (indent(), RG[0], RG[1], RG[2], RG[3]))

            print()

            print("%sMatrix Lmat(%s->extent(0), %s->extent(1), true);" % (indent(), lgammas[0][0], lgammas[0][0]))
            print("%sMatrix Rmat(%s->extent(0), %s->extent(1), true);" % (indent(), rgammas[0][0], rgammas[0][0]))

            print()

            abc = "%s" % left_indices[0]
            for i in range(1, len(oprod.left)):
                abc += " + %s%s" % (left_indices[i], ("*" + lnorb)*i)

            pqr = "%s" % right_indices[0]
            for i in range(1,len(oprod.right)):
                pqr += " + %s%s" %(right_indices[i], ("*" + rnorb)*i)

            integral_strings = [x.integrals(mo2e, input_indices[0:ninput], left_indices[0:len(x.left)], right_indices[0:len(x.right)]) for x in oprodvec]

            if len(oprod.right) > len(oprod.left): # inner loop is over right terms
                # start loops over left
                for i in range(len(oprod.left)):
                    print("%sfor (int %s = 0; %s < %s; ++%s) {" % (indent(), left_indices[i], left_indices[i], lnorb, left_indices[i]))
                    open_code_block()
                print("%sLmat.zero();" % indent())
                print("%sblas::ax_plus_y_n(1.0, &(*Lgamma)(0, 0, %s), Lmat.size(), Lmat.data());" % (indent(), abc))
                print()

                print("%sRmat.zero();" % indent())

                # start loops over right
                for i in range(len(oprod.right)):
                    print("%sfor (int %s = 0; %s < %s; ++%s) {" % (indent(), right_indices[i], right_indices[i], rnorb, right_indices[i]))
                    open_code_block()

                for i, RG in enumerate(rgammas):
                    print("%sblas::ax_plus_y_n(%1.1f * %s, &(*%s)(0, 0, %s), Rmat.size(), Rmat.data());" % (indent(), oprodvec[i].factor(), integral_strings[i], RG[0], pqr))

                # close right loops
                for i in range(len(oprod.right)):
                    close_code_block()

                rtrans = "true" if oprod.right.conj else "false"
                ltrans = "true" if oprod.left.conj else "false"

                phase = "left_phase" if len(oprod.right)%2==1 else 1.0

                print("%skronecker_product(%s, %s, Rmat, %s, Lmat, *out_block);" % (indent(), phase, rtrans, ltrans))

                for i in range(len(oprod.left)):
                    close_code_block()

            else: #inner loop is over left terms
                # write loops over Right terms
                for i in range(len(oprod.right)):
                    print("%sfor (int %s = 0; %s < %s; ++%s) {" % (indent(), right_indices[i], right_indices[i], rnorb, right_indices[i]))
                    open_code_block()

                print("%sRmat.zero();" % indent())
                print("%sblas::ax_plus_y_n(1.0, &(*Rgamma)(0, 0, %s), Rmat.size(), Rmat.data());" % (indent(), pqr))
                print()
                print("%sLmat.zero();" % indent())
                # loops over Left terms
                for i in range(len(oprod.left)):
                    print("%sfor (int %s = 0; %s < %s; ++%s) {" % (indent(), left_indices[i], left_indices[i], lnorb, left_indices[i]))
                    open_code_block()

                for i, LG in enumerate(lgammas):
                    print("%sblas::ax_plus_y_n(%1.1f * %s, &(*%s)(0, 0, %s), Lmat.size(), Lmat.data());" % (indent(), oprodvec[i].factor(), integral_strings[i], LG[0], abc))

                # close Left loops
                for i in range(len(oprod.left)):
                    close_code_block()

                rtrans = "true" if oprod.right.conj else "false"
                ltrans = "true" if oprod.left.conj else "false"

                phase = "left_phase" if len(oprod.right)%2==1 else 1.0

                print("%skronecker_product(%s, %s, Rmat, %s, Lmat, *out_block);" % (indent(), phase, rtrans, ltrans))

                # close Right loops
                for i in range(len(oprod.right)):
                    close_code_block()

            if nops > 1:
                close_code_block()
                #print()

        print("\n%s// add to map if large enough" % indent())
        print("%sif (out_block->rms() > thresh_)" % indent())
        print("%s  out.emplace(make_pair<size_t, size_t>(%s.offset, spair.offset), out_block);" % (indent(), tpair))

        if not both_diag:
            close_code_block()

        close_code_block()
        #print()

    close_code_block()
    print()
    print("%sreturn make_shared<BlockSparseMatrix>(%s.nstates, %s.nstates, out);" % (indent(), tinfo, sinfo))
    close_code_block()


generate_operator("S_a", [ContractedOperator([CreaA(), CreaA(), AnniA()], "S", 1.0), ContractedOperator([CreaA(), CreaB(), AnniB()], "S", 1.0)], 1)
print()
print()
generate_operator("S_b", [ContractedOperator([CreaB(), CreaB(), AnniB()], "S", 1.0), ContractedOperator([CreaB(), CreaA(), AnniA()], "S", 1.0)], 1)
print()
print()
generate_operator("Q_aa", [ContractedOperator([CreaA(), AnniA()], "QJK", 1.0), ContractedOperator([CreaB(), AnniB()], "QJ", 1.0)], 2)
print()
print()
generate_operator("Q_bb", [ContractedOperator([CreaB(), AnniB()], "QJK", 1.0), ContractedOperator([CreaA(), AnniA()], "QJ", 1.0)], 2)
print()
print()
generate_operator("Q_ab", [ContractedOperator([CreaB(), AnniA()], "QK", -1.0)], 2)
print()
print()
generate_operator("P_aa", [ContractedOperator([AnniA(), AnniA()], "P", 0.5)], 2)
print()
print()
generate_operator("P_bb", [ContractedOperator([AnniB(), AnniB()], "P", 0.5)], 2)
print()
print()
generate_operator("P_ab", [ContractedOperator([AnniB(), AnniA()], "P", 1.0)], 2)
