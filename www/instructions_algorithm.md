# Algorithm

For detailed description of the algorithm and its development see (Kovaleva et al 2023).

The algorithm searches for an arrangement of addresses <i>A</i>:</p>
<p>
  <math>
    <mi>A</mi>
    <mo>=</mo>
    <mo>{</mo>
    <msub><mi>a</mi><mn>1</mn></msub>,
    <msub><mi>a</mi><mn>2</mn></msub>,
    <msub><mi>a</mi><mn>3</mn></msub>,
    <mo>&#x2026;</mo>, <!-- Ellipsis -->
    <msub><mi>a</mi><mi>M</mi></msub>
    <mo>}</mo>
  </math>
<p>where <i>M</i> is the total number of tested peptides. Each address <i>a<sub>j</sub></i> corresponds to a peptide <i>E<sub>j</sub></i>.

The produced arrangement should meet all criteria described above.</p>

<p>Each address is encoded as a binary string with <i>n</i> bits, where <i>n</i> is equal to the total number of pools. Each digit <i>b<sub>i</sub></i> in an address <i>a<sub>j</sub></i> corresponds to a pool <i>p<sub>i</sub></i>, where <i>b<sub>i</sub> = 1</i> means the presence of a peptide <i>E<sub>j</sub></i> in pool <i>p<sub>i</sub></i>, and <i>b<sub>i</sub> = 0</i> means the absence of a peptide <i>E<sub>j</sub></i> in pool <i>p<sub>i</sub></i>. Two consecutive addresses <i>a<sub>j</sub></i> and <i>a<sub>j+1</sub></i> form a union <i>u<sub>j</sub></i>:</p>
<math>
  <mi>a</mi><sub><mi>j</mi></sub>
  <mo>&#x222A;</mo> <!-- Union symbol -->
  <mi>a</mi><sub><mi>j+1</mi></sub>
  <mo>=</mo>
  <mi>u</mi><sub><mi>j</mi></sub>
</math>
</p>
<p>Together, alternating addresses and unions form an extended arrangement <i>AU</i>:</p>
<math>
  <mi>AU</mi>
  <mo>=</mo>
  <mo>{</mo>
  <mi>a</mi><sub><mn>1</mn></sub>, <mi>u</mi><sub><mn>1</mn></sub>, <mi>a</mi><sub><mn>2</mn></sub>, <mi>u</mi><sub><mn>2</mn></sub>, <mi>a</mi><sub><mn>3</mn></sub>, <mi>u</mi><sub><mn>3</mn></sub>, <mo>&#x2026;</mo>, <!-- Ellipsis -->
  <mi>a</mi><sub><mi>M-1</mi></sub>, <mi>u</mi><sub><mi>M-1</mi></sub>, <mi>a</mi><sub><mi>M</mi></sub>
  <mo>}</mo>
</math>
</p>
<p>The union <i>u<sub>j</sub></i> also is encoded as a binary string with each digit <i>b<sub>i</sub></i> corresponding to a pool <i>p<sub>i</sub></i>. <i>b<sub>i</sub> = 1</i> in <i>u<sub>j</sub></i> means that <i>overlap(E<sub>j</sub>, E<sub>j+1</sub>)</i> is present in pool <i>p<sub>i</sub></i>, and <i>b<sub>i</sub> = 0</i> indicates its absence in pool <i>p<sub>i</sub></i>.</p>

<p>The arrangement of addresses <i>A</i> is characterized by its balance:</p>
<math>
  <mi>B</mi>
  <mo>=</mo>
  <mo>{</mo>
  <mi>w</mi><sub><mn>1</mn></sub>,
  <mi>w</mi><sub><mn>2</mn></sub>,
  <mi>w</mi><sub><mn>3</mn></sub>,
  <mo>&#x2026;</mo>, <!-- Ellipsis -->
  <mi>w</mi><sub><mi>n</mi></sub>
  <mo>}</mo>
</math>

<p>where <i>n</i> is the total number of pools, and <i>w<sub>i</sub></i> represents the total number of peptides in pool <i>p<sub>i</sub></i> and is calculated as:</p>

<math>
  <mi>w</mi>
  <msub>
  	<none/>
    <mi>i</mi>
  </msub>
  <mo>=</mo>
  <munderover>
    <mo>&#x2211;</mo> <!-- Summation symbol -->
    <mrow>
      <mi>k</mi>
      <mo>=</mo>
      <mn>1</mn>
    </mrow>
    <mi>M</mi>
  </munderover>
  <msub>
    <mi>a</mi>
    <mi>k</mi>
  </msub>
  <mo>[</mo>
  <mi>i</mi>
  <mo>]</mo>
</math>

<p>where <i>M</i> is the total number of peptides, <i>a<sub>k</sub></i> is an address, and <i>a<sub>k</sub>[i]</i> is its <i>i</i>th bit.</p>

The algorithm works as follows:

- initiate with arbitrary <i>a<sub>1</sub></i> and <i>u<sub>1</sub></i>:<p>
<math>
  <mrow>
    <mi>A</mi>
    <mi>U</mi>
    <mo>=</mo>
  </mrow>
  <mrow>
    <mo form="prefix" stretchy="false">{</mo>
    <msub>
      <mi>a</mi>
      <mn>1</mn>
    </msub>
    <mo separator="true">,</mo>
    <msub>
      <mi>u</mi>
      <mn>1</mn>
    </msub>
    <mo form="postfix" stretchy="false">}</mo>
  </mrow>
</math>

- produce all possible next addresses that can together with <i>a<sub>1</sub></i> form <i>u<sub>1</sub></i>, and discard addresses already present in <i>A</i>;

- calculate the distortion <i>D</i> to balance <i>B</i> introduced by adding each address to <i>A</i> with the addition of a small random number <math>
  <mi>ϵ</mi>
</math>:<p>
	<math>
  <mrow>
    <mi>D</mi>
    <mo>=</mo>
  </mrow>
  <mrow>
    <mtext>Var</mtext>
    <mo form="prefix" stretchy="false">(</mo>
    <msub>
      <mi>B</mi>
      <mrow>
        <mi>n</mi>
        <mi>e</mi>
        <mi>w</mi>
      </mrow>
    </msub>
    <mo form="postfix" stretchy="false">)</mo>
    <mo>−</mo>
  </mrow>
  <mrow>
    <mtext>Var</mtext>
    <mo form="prefix" stretchy="false">(</mo>
    <mi>B</mi>
    <mo form="postfix" stretchy="false">)</mo>
    <mo>+</mo>
  </mrow>
  <mrow>
    <mi>ϵ</mi>
  </mrow>
	</math><p>
	where <math>
  <mrow>
    <mn>0.0</mn>
    <mo>≤</mo>
  </mrow>
  <mrow>
    <mi>ϵ</mi>
    <mo>&lt;</mo>
  </mrow>
  <mrow>
    <mn>1.0</mn>
  </mrow>
</math>

- sort the list of possible addresses based on <i>D</i> and add the first one to the arrangement <i>A</i> and to the extended arrangement <i>AU</i>:<p>
	<math>
  <mrow>
    <mi>A</mi>
    <mo>=</mo>
  </mrow>
  <mrow>
    <mo form="prefix" stretchy="false">{</mo>
    <msub>
      <mi>a</mi>
      <mn>1</mn>
    </msub>
    <mo separator="true">,</mo>
    <msub>
      <mi>a</mi>
      <mn>2</mn>
    </msub>
    <mo form="postfix" stretchy="false">}</mo>
  </mrow>
</math><p>
	<math>
  <mrow>
    <mi>A</mi>
    <mi>U</mi>
    <mo>=</mo>
  </mrow>
  <mrow>
    <mo form="prefix" stretchy="false">{</mo>
    <msub>
      <mi>a</mi>
      <mn>1</mn>
    </msub>
    <mo separator="true">,</mo>
    <msub>
      <mi>u</mi>
      <mn>1</mn>
    </msub>
    <mo separator="true">,</mo>
    <msub>
      <mi>a</mi>
      <mn>2</mn>
    </msub>
    <mo form="postfix" stretchy="false">}</mo>
  </mrow>
</math>

- take the two last components of <i>AU</i> (<math>
  <msub>
    <mi>u</mi>
    <mn>1</mn>
  </msub>
</math> and <math>
  <msub>
    <mi>a</mi>
    <mn>2</mn>
  </msub>
</math>), generate  a list of all possible unions based thereon and discard unions already present in <i>AU</i>;

- calculate the distortion <i>D</i> to balance <i>B</i> introduced by adding each union to <i>A</i> with the addition of a small random number <math>
  <mi>ϵ</mi>
</math>;

- sort the list of possible unions based on <i>D</i> and add the first one to the extended arrangement <i>AU</i>:<p>
<math>
  <mrow>
    <mi>A</mi>
    <mi>U</mi>
    <mo>=</mo>
  </mrow>
  <mrow>
    <mo form="prefix" stretchy="false">{</mo>
    <msub>
      <mi>a</mi>
      <mn>1</mn>
    </msub>
    <mo separator="true">,</mo>
    <msub>
      <mi>u</mi>
      <mn>1</mn>
    </msub>
    <mo separator="true">,</mo>
    <msub>
      <mi>a</mi>
      <mn>2</mn>
    </msub>
    <mo separator="true">,</mo>
    <msub>
      <mi>u</mi>
      <mn>2</mn>
    </msub>
    <mo form="postfix" stretchy="false">}</mo>
  </mrow>
</math>

- take the two last components of <i>AU</i> and continue with 2., which leads to traversing through space of all possible addresses and unions;

- stop when the length of the arrangement <i>A</i> is equal to <i>M</i>.

In copepodTCR, the algorithm is implemented recursively, allowing for backtracking if it reaches a dead end, enabling it to select another address in the list of next addresses, and another union in the list of next unions.
