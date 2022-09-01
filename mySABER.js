// Matrix creation of size x*y
function createMatrix(x, y)
{
    let matrix =  new Array(x);
    for(let i = 0; i < x; i++)
    {
        matrix[i] = new Array(y);
    }
    return matrix;
}

// Random Matrix creation, each element is chosen at random
function createRandomMatrix(x, y, n, z)
{
    let matrix = new Array(x);
    for(let i = 0; i < x; i++)
    {
        matrix[i] =  new Array(y);
        for(let j = 0; j < y; j++)
        {
            matrix[i][j] = new Array(n);
            for(let r = 0; r < n; r++)
            {
                matrix[i][j][r] = nextInt(z);
            }
        }
    }
    return matrix;
}

// Generating the secret matrix s (INCOMPLETE)
function genSecret()
{
    let secret = new Array(k);
    let fours = new Array(n);
    fours.fill(4); 

    for(let i = 0; i < k; i++)
    {
        secret[i] = new Array(n);
    }

    var binom = SJS.Binomial(8, 0.5);

    for(let i = 0; i < k; i++)
    {
        secret[i] = substractArrays(binom.sample(n), fours);
    }

    return secret;
}

// This function returns the Hamming weight of the input bit string aka the number of nonzero elements
function HammightWeight(bitstring)
{
    let counter = 0;

    for(let i = 0; i < bitstring.length; i++)
    {
        if(bitstring[i] != 0)
        {
            counter++;
        }
    }

    return counter;
}

// Functions on Matrixes:
// Returns the transpose of a matrix A
function transpose(matrix) 
{
    output = matrix[0].map((_, colIndex) => matrix.map(row => row[colIndex]));
    return output;
}

// Returns the matrix A modulo q
function mod(A, q) 
{
	let x = A.length;
	let y = A[0].length;
    let B = JSON.parse(JSON.stringify(A));

	if (q <= 0) 
    {
		alert("Modulus not positive");
		return;
	}
	for (let i = 0; i < x; i++) 
    {
		for (let j = 0; j < y; j++) 
        {
			B[i][j] = ((A[i][j] % q) + q) % q;
		}
	}

    return B;
} 

// Returns the vector v modulo q
function modvector(v, q)
{
    let x = v.length;
    let vector = JSON.parse(JSON.stringify(v));

    for(let i = 0; i < x; i++)
    {
        vector[i] = ((v[i] % q) + q) % q
    }

    return vector;
}

// Multiplies vector 
function multiplyVectors(a, b)
{
    let result = new Array(a[0].length).fill(0);
    
    for(let i = 0; i < a.length; i++)
    {
        for(let j = 0; j < a[0].length; j++)
        {
            for(let l = 0; l < b[0].length; l++)
            {
                result[j+l] += a[i][j]*b[i][l];
            }
        }
    }

    return result;
}

// Returns the right side multiplication of a Matrix A with a vector v, M = A*v
function multiplyMatrixVector(A, v)
{
    let x = A.length;
    let y = A[0].length;
    let z = v.length;

    if(y != z)
    {
        alert("Inner dimensions have to be the same");
        return;
    }

    let M = new Array(z);

    for(let i = 0; i < x; i++)
    {
        M[i] = 0;
        for(let j = 0; j < y; j++)
        {
            M[i] += A[i][j]*v[j];
        }
    }
    return M;
}

// Create constant matrix
function createConstantVector()
{
    let v = new Array(l);

    for(let i = 0; i < l; i++)
    { 
        v[i] = new Array(n);
        for(let j = 0; j < n; j++)
        {
            v[i][j] = 2 ** (eq - ep - 1);
        }
    }

    return v;
}

// Creating the constant vector h2
function createConstantVector2()
{
    let v = new Array(n);

    for(let i = 0; i < n; i++)
    {
        v[i] = 2 ** (ep - 2) - 2 ** (ep - et - 1) + 2 ** (eq - ep - 1);
    }

    return v;
}

// Multiplication of A and s
function MatrixVectorMul(M, v, q)
{
    let result = new Array(k);
    let c = new Array(n);

    for(let i = 0; i < k; i++)
    {
        result[i] = new Array(n);
        result[i].fill(0);
    }

    for(let i = 0; i < l; i++)
    {
        c.fill(0);

        for(let j = 0; j < l; j++)
        {
            c = additionArrays(c, PolyMul(M[i][j], v[j], q));
        }
        result[i] = modvector(c, q);
    }

    return result;
}

// Matrix Addition v = a + b mod q
function additionVectorsmod(a, b, q)
{
    let v = new Array(a.length);

    if(a.length != b.length)
    {
        alert("Vector sizes not the same");
        return;
    }
    for(let i = 0; i < a.length; i++)
    {
        v[i] = new Array(n);

        for(let j = 0; j < n; j++)
        {
            v[i][j] = (((a[i][j] + b[i][j]) % q) + q)%q;
        }
    }
    
    return v;
}

// Adding array a to array b
function additionArrays(a, b)
{
    let x = a.length;
    let y = b.length;
    let c = new Array(a.length);

    if(x != y)
    {
        alert("Array size missmatch");
    }
    for(let i = 0; i < a.length; i++)
    {
        c[i] =  a[i] + b[i];
    }

    return c;
}

// Substract array b from array a
function substractArrays(a, b)
{
    let x = a.length;
    let y = b.length;
    let c = new Array(a.length);

    if(x != y)
    {
        alert("Array size missmatch");
    }

    for(let i = 0; i < a.length; i++)
    {
        c[i] = a[i] - b[i];
    }
    return c;
}

// Takes the coefficients of the polynomial 'a' and shifts them right by q - p positions
function shiftRight(a, eq)
{
    let v = new Array(a.length);

    for(let i = 0; i < a.length; i++)
    {
        v[i] = a[i] >> (eq);
    }

    return v;
}

// Takes the coefficients of the polynomial 'a' and shifts them left by q - p positions
function shiftLeft(a, eq)
{
    let v = new Array(a.length);

    for(let i = 0; i < a.length; i++)
    {        
        v[i] = a[i] << eq;
    }

    return v;
}

// This function performs polynomial multiplication for polynomials a and b modulo p with product c
function PolyMul(a, b, p)
{
    let c = new Array(n);
    c.fill(0);

    for(let i = 0; i < n; i++)
    {
        for(let j = 0; j < n; j++)
        {
            c[(i + j)%n] += a[i]*b[j];
        }
    }
    for(let i = 0; i < c.length; i++)
    {
        c[i] = ((c[i] % p) + p) % p;
    }

    return c;
}

// This function takes a vector v1 and a vector v2 and computes the inner product modulo p which is a polynomial c
function Innerprod(v1, v2, p)
{
    let c = new Array(n);
    c.fill(0);

    for(let i = 0; i < l; i++)
    {
        c = additionArrays(c, PolyMul(v1[i], v2[i], p));
    }
    for(let i = 0; i < c.length; i++)
    {
        c[i] = ((c[i] % p) + p) % p;
    }

    return c;
}

// This functions takes a polynomial from RN and transforms it into a byte string of length k x 256 / 8 
function pol2BS(pol, k)
{
    let bits = new Array(pol.length);
    let conc = new Array(pol.length + 1);

    for(let i = 0; i < pol.length; i++)
    {
        bits[i] = new Array(k);

        for(let j = 0; j < k; j++)
        {
            bits[i][k - j - 1] = ( pol[i] >> j ) & 1; // Create the bits from the coefficients of the polynomial
        }
        conc[i] = bits[i].join(""); // Create the concatenated bit string
    }

    let bitstring = conc.join(""); 

    /* Turn the bit string into a byte string (every 8 bits become a byte aka every 4 bits become one hexadecimal)
    let bytestring = new Array(bitstring.length/4);
    let counter = 0;

    for(let i = 0; i < bitstring.length; i+= 8)
    {
        byte[0] = bitstring[i];
        byte[1] = bitstring[i + 1];
        byte[2] = bitstring[i + 2];
        byte[3] = bitstring[i + 3];
        byte[4] = bitstring[i + 4];
        byte[5] = bitstring[i + 5];
        byte[6] = bitstring[i + 6];
        byte[7] = bitstring[i + 7];
        let firstHex = byte[0] + "" + byte[1] + "" + byte[2] + "" + byte[3];
        let secondHex = byte[4] + "" + byte[5] + "" + byte[6] + "" + byte[7];
        firstHex = parseInt(firstHex, 2).toString(16);
        secondHex = parseInt(secondHex, 2).toString(16);
        bytestring[counter] = firstHex;
        counter++;
        bytestring[counter] = secondHex;
        counter++;
    } */

    return bitstring;
}

// This function takes a vector from R^lx1 and transforms it into a byte string of length l x k x 256 / 8 
function polvec2BS(v, p)
{
    let bytestrings =  new Array(v.length);

    for(i = 0; i < v.length; i++)
    {
        bytestrings[i] =  pol2BS(v[i], p);
    }

    return bytestrings.join("");
}

// This function takes a byte string of length k x 256 / 8 and transforms it into a polynomial 
function bs2pol(bitstring, k)
{
    let chunks = new Array(k);
    let coeff = new Array(bitstring.length/k);

    for(let i = 0; i < bitstring.length/k; i += 1)
    {
        for(let j = 0; j < k; j++)
        {
            chunks[j] = bitstring[i*k + j];
        }
        coeff[i] = parseInt(chunks.join(""), 2).toString(10);
    }

    /* for(let i = 0; i < n*k; i += 3)
    {
        coeff[i/3] = parseInt((bitstring[i] + "" + bitstring[i + 1] + "" + bitstring[i + 2]), 2).toString(2);
        coeff[i/3] = coeff[i/3].padStart(3, '0');
    } */ 

    return coeff;
}

// This function takes a byte string of length l x k x 256/8 and transforms it into a vector
function bs2polvec(bitstring, p)
{
    let v = new Array(l);
    let chunk = bitstring.length/3; 
    
    for(let i = 0; i < l; i++)
    {
        v[i] = bs2pol(bitstring.slice(i*chunk, (i + 1)*chunk), p);
    }   

    return v;
}

// This function generates a random 256 bit message we use for encryption and decryption
function genMessage(size)
{
    let bitstring = new Array(size);

    for(let i = 0; i < size; i++)
    {
        bitstring[i] = Math.floor(Math.random() * 2);
    }

    /* for(let i = 0; i < bitstring.length*k; i+= 8) HEXADECIMALS
    {
        byte[0] = bitstring[i];
        byte[1] = bitstring[i + 1];
        byte[2] = bitstring[i + 2];
        byte[3] = bitstring[i + 3];
        byte[4] = bitstring[i + 4];
        byte[5] = bitstring[i + 5];
        byte[6] = bitstring[i + 6];
        byte[7] = bitstring[i + 7];
        let firstHex = byte[0] + "" + byte[1] + "" + byte[2] + "" + byte[3];
        let secondHex = byte[4] + "" + byte[5] + "" + byte[6] + "" + byte[7];
        firstHex = parseInt(firstHex, 2).toString(16);
        secondHex = parseInt(secondHex, 2).toString(16);
        bytestring[counter] = firstHex;
        counter++;
        bytestring[counter] = secondHex;
        counter++; 
    } */

    return bitstring.join('');
}

// Verifies if the message at the end of the process is equal to the message at the start
function verify(a, b)
{
    for(let i = 0; i < n; i++)
    {
        if(a[i] !== b[i])
        {
            console.log("The messages are not the same");
            return false;
        }
    }

    console.log("The messages are the same!");
    return true;
}

// ------------------------KEY-ENCAPSULATION-MECHANISM--------------------------------
// Key Generation function

function keyGen()
{
    console.log("KEY GENERATION--------------------------");
    console.log("Creating the PublicKey and SecretKey");
    A = createRandomMatrix(k, k, n, q); // k*k mod q
    var h = createConstantVector(); // k*1 mod q
    var b = new Array(k);
    var bp = new Array(k);
    let s = genSecret();

    for(let i = 0; i < k; i++)
    {
        b[i] = new Array(n);
    }

    b = additionVectorsmod(MatrixVectorMul(transpose(A) , s, q), h, q); 

    for(let i = 0; i < l; i++)
    {
        bp[i] = shiftRight(b[i], eq - ep);
    } 

    // Public key pk and secret key sk
    sk = polvec2BS(s, eq); // This is a bitstring with size 8 x SABER_INDCPA_SECRETKEYBYTES = 9984
    pk = polvec2BS(bp, ep); // This is a bitstring with size ( SABER_INDCPA_SECRETKEYBYTES - SABER_SEEBYTES ) x 8 = 7680
}

// Encryption
function encrypt()
{
    console.log("ENCRYPTION------------------------------");
    console.log("Input: Message m, seed(s2) and PublicKey ");
    var s2 = genSecret();
    message = genMessage(n); //  message m is 256 bits
    var mp =  new Array(n); // polynomial with bytes as coefficients containing the message

    let b2 = additionVectorsmod(MatrixVectorMul(A, s2, q), h, q); // A'*s + h

    for(let i = 0; i < l; i++)
    {
        b2[i] = shiftRight(b2[i], eq - ep); // shift right by (eq - ep)
    }

    b = bs2polvec(pk, ep);
    v2 = Innerprod(b, s2, p);
    mp = bs2pol(message, 1);
    mp = shiftLeft(mp, ep - 1);
    cm = shiftRight(modvector(additionArrays(substractArrays(v2, mp) , h1), p), ep - et);
    ciphertext =  pol2BS(cm, et).concat(polvec2BS(b2, ep));
}

//Decryption
function decrypt()
{
    let b2, v, cm, ct;
    console.log("DECRYPTION------------------------------");
    console.log("Input: CipherText and SecretKey");
    s = bs2polvec(sk, eq); // Turn secret key into a polynomial vector s

    // ( cm || ct) =  CipherText
    cm = ciphertext.slice(0, n*et);
    ct = ciphertext.slice(n*et);

    cm = bs2pol(cm, et);
    cm = shiftLeft(cm, ep - et);
    b2 = bs2polvec(ct, ep);
    v = Innerprod(b2, mod(s, p), p);

    message2 = shiftRight(modvector(additionArrays(substractArrays(v, cm), h2), p)  , ep - 1);
    message2 = pol2BS(message2, 1);

    console.log("ORIGINAL MESSAGE: " + message);
    console.log("CREATED MESSAGE: " + message2);
}

//Returns the next pseudorandom, uniformly distributed integer between 0(inclusive) and q-1(inclusive)
function nextInt(q) 
{
	return Math.floor(random() * q);
}

//Returns the pseudorandom integer value between low(inclusive) and high(inclusive)
function rangeValue(low, high) 
{
	return Math.floor(random() * (high - low + 1) + low);
}

//-----------------------START-SABER--------------------------
// defining variables
var q = 8192, p = 1024, k = 3, t= 16, m = 8, n = 256, l = 3; // Our base variables for dimensions
var eq = Math.log(q)/Math.log(2); ep = Math.log(p)/Math.log(2), et = Math.log(t)/Math.log(2); // The powers of our base variables
var pk, sk; // Public and Secret key
var A; // Our generated matrix
var ciphertext; // The ciphertext created in the encryption process
var h = createConstantVector();
var h1 = createConstantVector()[0]; // First element from h is the constant polynomial h1
var h2 = createConstantVector2(); // Creating the 2nd constant vector used for rounding operations
var message, message2; // Message used in encryption and message2 used in decryption

function testSaber()
{
    console.log("Testing SABER KEM ");
    console.log("Our variables are: ");
    console.log("q = " + q);
    console.log("p = " + p);
    console.log("k = " + k);
    console.log("t = " + t);
    console.log("μ = " + m);
    console.log("n = " + n);
    console.log("l = " + l);
    keyGen();
    encrypt();
    decrypt();

    verify(message, message2);
}

testSaber();