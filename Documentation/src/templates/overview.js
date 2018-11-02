module.exports = `
<h3 id="overview">Overview</h3>
<p>Optimizer provides C++ solutions to optimization problems.  The goal is to build an Optimization library for C++ which can provide solutions for single variable and multi-variable optimization, constrained and unconstrained problems. The users should be able to select the algorithm to be used and also access any intermediate data which is obtained whilst running the algorithm.</p>
<h3 id="a-simple-first-program">A simple first program</h3>
<pre><code><span class="hljs-meta">#<span class="hljs-meta-keyword">include</span> <span class="hljs-meta-string">&lt;iostream&gt;</span></span>
<span class="hljs-meta">#<span class="hljs-meta-keyword">include</span> <span class="hljs-meta-string">&lt;Eigen/Dense&gt;</span></span>
<span class="hljs-meta">#<span class="hljs-meta-keyword">include</span> <span class="hljs-meta-string">"../Optimizer/optimizer"</span></span>
<span class="hljs-keyword">using</span> <span class="hljs-keyword">namespace</span> <span class="hljs-built_in">std</span>;
<span class="hljs-keyword">using</span> <span class="hljs-keyword">namespace</span> Eigen;
<span class="hljs-keyword">using</span> <span class="hljs-keyword">namespace</span> Optimizer;

<span class="hljs-function"><span class="hljs-keyword">int</span> <span class="hljs-title">main</span> <span class="hljs-params">()</span> </span>{
  <span class="hljs-built_in">cout</span> &lt;&lt; <span class="hljs-string">"Using Function: (x + 10)^2 for single variable algorithms testing."</span> &lt;&lt; <span class="hljs-built_in">endl</span>;
  <span class="hljs-keyword">double</span> ipt = <span class="hljs-number">5.4</span>;

  <span class="hljs-built_in">cout</span> &lt;&lt; <span class="hljs-string">"Test Bounding Phase:"</span> &lt;&lt; <span class="hljs-built_in">endl</span>;
  Vector2d range = BoundingPhase(func, ipt);
  <span class="hljs-built_in">cout</span> &lt;&lt; <span class="hljs-string">"Range from bounding Phase for initial point :"</span> &lt;&lt; ipt &lt;&lt; <span class="hljs-built_in">endl</span>;
  <span class="hljs-built_in">cout</span> &lt;&lt; range &lt;&lt; <span class="hljs-built_in">endl</span>;

  <span class="hljs-built_in">cout</span> &lt;&lt; <span class="hljs-string">"Test Exhaustive Search:"</span> &lt;&lt; <span class="hljs-built_in">endl</span>;
  range = Exhaustive(func, ipt);
  <span class="hljs-built_in">cout</span> &lt;&lt; <span class="hljs-string">"Range from Exhaustive Search for initial point :"</span> &lt;&lt; ipt &lt;&lt; <span class="hljs-built_in">endl</span>;
  <span class="hljs-built_in">cout</span> &lt;&lt; range &lt;&lt; <span class="hljs-built_in">endl</span>;

  <span class="hljs-built_in">cout</span> &lt;&lt; <span class="hljs-string">"Derivatives at "</span> &lt;&lt; ipt &lt;&lt; <span class="hljs-built_in">endl</span>;
  <span class="hljs-built_in">cout</span> &lt;&lt; Derivative(func, ipt) &lt;&lt; <span class="hljs-built_in">endl</span>;

  <span class="hljs-built_in">cout</span> &lt;&lt; <span class="hljs-string">"Finding optimal point using above range for Newton Rapshon Method."</span> &lt;&lt; <span class="hljs-built_in">endl</span>;
  <span class="hljs-built_in">cout</span> &lt;&lt; <span class="hljs-string">"Optimal Point is: "</span>;
  <span class="hljs-built_in">cout</span> &lt;&lt; NewtonRapshon (func, range) &lt;&lt; <span class="hljs-built_in">endl</span>;
}</code></pre>
`