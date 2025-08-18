<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.12.0">
  <compound kind="file">
    <name>convert_to_layered_sparse.hpp</name>
    <path>tatami_layered/</path>
    <filename>convert__to__layered__sparse_8hpp.html</filename>
    <class kind="struct">tatami_layered::ConvertToLayeredSparseOptions</class>
    <namespace>tatami_layered</namespace>
  </compound>
  <compound kind="file">
    <name>read_layered_sparse_from_matrix_market.hpp</name>
    <path>tatami_layered/</path>
    <filename>read__layered__sparse__from__matrix__market_8hpp.html</filename>
    <class kind="struct">tatami_layered::ReadLayeredSparseFromMatrixMarketOptions</class>
    <namespace>tatami_layered</namespace>
  </compound>
  <compound kind="file">
    <name>tatami_layered.hpp</name>
    <path>tatami_layered/</path>
    <filename>tatami__layered_8hpp.html</filename>
    <namespace>tatami_layered</namespace>
  </compound>
  <compound kind="struct">
    <name>tatami_layered::ConvertToLayeredSparseOptions</name>
    <filename>structtatami__layered_1_1ConvertToLayeredSparseOptions.html</filename>
    <member kind="variable">
      <type>std::size_t</type>
      <name>chunk_size</name>
      <anchorfile>structtatami__layered_1_1ConvertToLayeredSparseOptions.html</anchorfile>
      <anchor>afac100e950aad3b286792c2a614c6cb7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_threads</name>
      <anchorfile>structtatami__layered_1_1ConvertToLayeredSparseOptions.html</anchorfile>
      <anchor>af699bfb9473eeae87cceade395a47d4b</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami_layered::ReadLayeredSparseFromMatrixMarketOptions</name>
    <filename>structtatami__layered_1_1ReadLayeredSparseFromMatrixMarketOptions.html</filename>
    <member kind="variable">
      <type>std::size_t</type>
      <name>chunk_size</name>
      <anchorfile>structtatami__layered_1_1ReadLayeredSparseFromMatrixMarketOptions.html</anchorfile>
      <anchor>a5c30c7f13ca4d3886180bf0aa2fef730</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::size_t</type>
      <name>buffer_size</name>
      <anchorfile>structtatami__layered_1_1ReadLayeredSparseFromMatrixMarketOptions.html</anchorfile>
      <anchor>a157622423b2f19b451fa7700ab364e07</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>num_threads</name>
      <anchorfile>structtatami__layered_1_1ReadLayeredSparseFromMatrixMarketOptions.html</anchorfile>
      <anchor>af4e79dbd84ee6e472a22a3fd10fea03e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>tatami_layered</name>
    <filename>namespacetatami__layered.html</filename>
    <class kind="struct">tatami_layered::ConvertToLayeredSparseOptions</class>
    <class kind="struct">tatami_layered::ReadLayeredSparseFromMatrixMarketOptions</class>
    <member kind="function">
      <type>std::shared_ptr&lt; tatami::Matrix&lt; ValueOut_, IndexOut_ &gt; &gt;</type>
      <name>convert_to_layered_sparse</name>
      <anchorfile>namespacetatami__layered.html</anchorfile>
      <anchor>ae6f881696ce1e2fa7b9e2aeb3e38678b</anchor>
      <arglist>(const tatami::Matrix&lt; ValueIn_, IndexIn_ &gt; &amp;mat, const ConvertToLayeredSparseOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; tatami::Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>read_layered_sparse_from_matrix_market_text_file</name>
      <anchorfile>namespacetatami__layered.html</anchorfile>
      <anchor>a99e2feeee42e4ab42b8f7db3e0144419</anchor>
      <arglist>(const char *filepath, const ReadLayeredSparseFromMatrixMarketOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; tatami::Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>read_layered_sparse_from_matrix_market_some_file</name>
      <anchorfile>namespacetatami__layered.html</anchorfile>
      <anchor>a1cd8b6c524c08a2672f01fa188086509</anchor>
      <arglist>(const char *filepath, const ReadLayeredSparseFromMatrixMarketOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; tatami::Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>read_layered_sparse_from_matrix_market_gzip_file</name>
      <anchorfile>namespacetatami__layered.html</anchorfile>
      <anchor>a3b7a376824c0e0159f3c12ea24b49210</anchor>
      <arglist>(const char *filepath, const ReadLayeredSparseFromMatrixMarketOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; tatami::Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>read_layered_sparse_from_matrix_market_text_buffer</name>
      <anchorfile>namespacetatami__layered.html</anchorfile>
      <anchor>a5514376d8fe31e185285ddaadd763081</anchor>
      <arglist>(const unsigned char *contents, std::size_t length, const ReadLayeredSparseFromMatrixMarketOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; tatami::Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>read_layered_sparse_from_matrix_market_some_buffer</name>
      <anchorfile>namespacetatami__layered.html</anchorfile>
      <anchor>a135938c24f4fe88a7b3590ffeed4c2c8</anchor>
      <arglist>(const unsigned char *contents, std::size_t length, const ReadLayeredSparseFromMatrixMarketOptions &amp;options)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; tatami::Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>read_layered_sparse_from_matrix_market_zlib_buffer</name>
      <anchorfile>namespacetatami__layered.html</anchorfile>
      <anchor>a8e9ed17125a431ef82be62b389346126</anchor>
      <arglist>(const unsigned char *contents, std::size_t length, const ReadLayeredSparseFromMatrixMarketOptions &amp;options)</arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>tatami helpers for creating layered matrices</title>
    <filename>index.html</filename>
    <docanchor file="index.html">md__2github_2workspace_2README</docanchor>
  </compound>
</tagfile>
