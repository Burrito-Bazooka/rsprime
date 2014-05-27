import unittest
import itertools

from rsprime import PFint, Polynomial, RSCoder

PF59int = PFint(59)

class TestPFint(unittest.TestCase):
    def test_arithmetic(self):
        a = PF59int(3)
        b = PF59int(9)

        self.assertEqual(a + b, 12)
        self.assertEqual(b + a, 12)
        self.assertEqual(3 + b, 12)
        self.assertEqual(a + 9, 12)

        self.assertEqual(a - b, 53)
        self.assertEqual(3 - b, 53)
        self.assertEqual(a - 9, 53)
        self.assertEqual(b - a, 6)
        self.assertEqual(9 - a, 6)
        self.assertEqual(b - 3, 6)

        self.assertEqual(a * b, 27)
        self.assertEqual(b * a, 27)
        self.assertEqual(3 * b, 27)
        self.assertEqual(a * 9, 27)

        self.assertEqual(b * b.inverse(), 1)
        self.assertEqual(b / b, 1)

        self.assertEqual(b / a, 3)
        self.assertEqual(9 / a, 3)
        self.assertEqual(b / 3, 3)

        self.assertRaises(Exception, lambda: b**a)
        self.assertEqual(b**3, 21)
        self.assertRaises(Exception, lambda: a**b)
        self.assertEqual(a**9, 36)

        self.assertEqual(b.inverse(), 46)
        self.assertEqual(b * 46, 1)

    def test_fermats_theorem(self):
        for x in range(1,59):
            self.assertEqual(PF59int(x)**58, 1)

class TestRSverify(unittest.TestCase):
    def setUp(self):
        self.coder = RSCoder(59,58,46)

    def test_one(self):
        """Tests a codeword without errors validates"""
        code = self.coder.encode("1Ah56Cfe4SXA", nostrip=True)

        self.assertTrue(self.coder.verify(code))

    def test_two(self):
        """Verifies that changing any single character will invalidate the
        codeword"""
        code = self.coder.encode("123456789abcdefghijkmnpqrstuvwxyzA", nostrip=True)

        for i, c in enumerate(code):
            # Change the value at position i and verify that the code is not
            # valid
            # Change it to a 0, unless it's already a 0
            if self.coder.mapper.decode(c) == 0:
                c = self.coder.mapper.encode(1)
            else:
                c = self.coder.mapper.encode(0)
            bad_code = code[:i] + c + code[i+1:]

            self.assertFalse(self.coder.verify(bad_code))

class TestRSdecoding(unittest.TestCase):
    def setUp(self):
        self.coder = RSCoder(59,58,46)
        self.string = "818878"

        codestr = self.coder.encode(self.string, nostrip=True)

        self.code = codestr

    def test_strip(self):
        """Tests that the nostrip feature works"""
        otherstr = self.string.rjust(46, "0")

        codestr = self.coder.encode(otherstr, nostrip=True)

        self.assertEqual(58, len(codestr))
        
        # Decode with default behavior: stripping of leading null bytes
        decode = self.coder.decode(codestr)
        decode2 = (self.coder.decode(codestr[:5] + "1" + codestr[6:]))

        self.assertEqual(self.string, decode)
        self.assertEqual(self.string, decode2)

        # Decode with nostrip
        decode = self.coder.decode(codestr, nostrip=True)
        decode2 = self.coder.decode(codestr[:5] + "1" + codestr[6:], nostrip=True)

        self.assertEqual(otherstr, decode)
        self.assertEqual(otherstr, decode2)

    def test_noerr(self):
        """Make sure a codeword with no errors decodes"""
        decode = self.coder.decode(self.code)
        self.assertEqual(self.string, decode)

    def test_oneerr(self):
        """Change just one byte and make sure it decodes"""
        for i, c in enumerate(self.code):
            newch = self.coder.mapper.encode( (self.coder.mapper.decode(c)+1) % 59)
            r = self.code[:i] + newch + self.code[i+1:]

            decode = self.coder.decode(r)

            self.assertEqual(self.string, decode)


    def test_6err(self):
        """Tests if 16 byte errors still decodes"""
        errors = [5, 6, 12, 13, 38, 40]
        r = list(self.coder.mapper.decode(x) for x in self.code)

        for e in errors:
            r[e] = (r[e] + 1)

        r = "".join(self.coder.mapper.encode(x) for x in r)
        decode = self.coder.decode(r)
        self.assertEqual(self.string, decode)

    def test_17err(self):
        """Kinda pointless, checks that 17 errors doesn't decode.
        Actually, this could still decode by coincidence on some inputs,
        so this test shouldn't be here at all."""
        errors = [5, 6, 12, 13, 22, 38, 40, 42]
        r = list(self.coder.mapper.decode(x) for x in self.code)

        for e in errors:
            r[e] = (r[e] + 50) % 256

        r = "".join(self.coder.mapper.encode(x) for x in r)
        decode = self.coder.decode(r)
        self.assertNotEqual(self.string, decode)


class TestPFPoly(unittest.TestCase):
    """Tests that the Polynomial class works when given PFint objects
    instead of regular integers
    """
    def test_add(self):
        one = Polynomial(map(PF59int,     (1,3,5,1)))
        two = Polynomial(map(PF59int, (5,3,58,1,6,8)))

        r = one + two

        self.assertEqual(r.coefficients,(PF59int(5),PF59int(3),PF59int(0),PF59int(4),PF59int(11),PF59int(9)))

    def test_sub(self):
        one = Polynomial(map(PF59int,     (8,3,5,1)))
        two = Polynomial(map(PF59int, (5,3,1,1,6,8)))
        r = one - two
        self.assertEqual(r.coefficients, (PF59int(54),PF59int(56),PF59int(7),PF59int(2),PF59int(58),PF59int(52)))

    def test_mul(self):
        one = Polynomial(map(PF59int,     (8,3,5,1)))
        two = Polynomial(map(PF59int, (5,3,1,1,6,8)))
        r = one * two
        self.assertEqual(r.coefficients, (40,39,42,31,0,29,55,46,8))

    def test_div(self):
        one = Polynomial(map(PF59int,     (1,58)))
        two = Polynomial(map(PF59int, (1,0,58)))
        q, r = divmod(two,one)
        self.assertEqual(q.coefficients, (PF59int(1),PF59int(1)))
        self.assertEqual(r.coefficients, (0,))

        # Make sure they multiply back out okay
        self.assertEqual(q*one + r, two)

    def test_div_scalar(self):
        """Tests division by a scalar"""
        numbers = map(PF59int, (5,20,50,10,34,58,0,48,33,25,4,5,2))
        scalar = PF59int(17)

        poly = Polynomial(numbers)
        scalarpoly = Polynomial(x0=scalar)

        self.assertEqual(
                (poly // scalarpoly).coefficients,
                tuple(map(lambda x: x / scalar, numbers))
                )

    def test_div_scalar2(self):
        """Test that dividing by a scalar is the same as multiplying by the
        scalar's inverse"""
        a = Polynomial(map(PF59int, (5,3,1,1,6,8)))

        scalar = PF59int(50)

        self.assertEqual(
                a * Polynomial(x0=scalar),
                a // Polynomial(x0=scalar.inverse())
                )



class TestPolynomial(unittest.TestCase):
    def test_add_1(self):
        one = Polynomial((2,4,7,3))
        two = Polynomial((5,2,4,2))

        r = one + two

        self.assertEqual(r.coefficients, (7, 6, 11, 5))

    def test_add_2(self):
        one = Polynomial((2,4,7,3,5,2))
        two = Polynomial((5,2,4,2))

        r = one + two

        self.assertEqual(r.coefficients, (2,4,12,5,9,4))

    def test_add_3(self):
        one = Polynomial((7,3,5,2))
        two = Polynomial((6,8,5,2,4,2))

        r = one + two

        self.assertEqual(r.coefficients, (6,8,12,5,9,4))

    def test_mul_1(self):
        one = Polynomial((2,4,7,3))
        two = Polynomial((5,2,4,2))

        r = one * two

        self.assertEqual(r.coefficients,
                (10,24,51,49,42,26,6))

    def test_div_1(self):
        one = Polynomial((1,4,0,3))
        two = Polynomial((1,0,1))

        q, r = divmod(one, two)
        self.assertEqual(q, one // two)
        self.assertEqual(r, one % two)

        self.assertEqual(q.coefficients, (1,4))
        self.assertEqual(r.coefficients, (-1,-1))

    def test_div_2(self):
        one = Polynomial((1,0,0,2,2,0,1,2,1))
        two = Polynomial((1,0,-1))

        q, r = divmod(one, two)
        self.assertEqual(q, one // two)
        self.assertEqual(r, one % two)

        self.assertEqual(q.coefficients, (1,0,1,2,3,2,4))
        self.assertEqual(r.coefficients, (4,5))

    def test_div_3(self):
        # 0 quotient
        one = Polynomial((1,0,-1))
        two = Polynomial((1,1,0,0,-1))

        q, r = divmod(one, two)
        self.assertEqual(q, one // two)
        self.assertEqual(r, one % two)

        self.assertEqual(q.coefficients, (0,))
        self.assertEqual(r.coefficients, (1,0,-1))

    def test_div_4(self):
	# no remander
        one = Polynomial((1,0,0,2,2,0,1,-2,-4))
        two = Polynomial((1,0,-1))

        q, r = divmod(one, two)
        self.assertEqual(q, one // two)
        self.assertEqual(r, one % two)

        self.assertEqual(q.coefficients, (1,0,1,2,3,2,4))
        self.assertEqual(r.coefficients, (0,))

    def test_getcoeff(self):
        p = Polynomial((9,3,3,2,2,3,1,-2,-4))
        self.assertEqual(p.get_coefficient(0), -4)
        self.assertEqual(p.get_coefficient(2), 1)
        self.assertEqual(p.get_coefficient(8), 9)
        self.assertEqual(p.get_coefficient(9), 0)

if __name__ == "__main__":
    unittest.main()
