<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.9.8">
  <compound kind="file">
    <name>convert_to_layered_sparse.hpp</name>
    <path>tatami_layered/</path>
    <filename>convert__to__layered__sparse_8hpp.html</filename>
    <namespace>tatami_layered</namespace>
  </compound>
  <compound kind="file">
    <name>read_layered_sparse_from_matrix_market.hpp</name>
    <path>tatami_layered/</path>
    <filename>read__layered__sparse__from__matrix__market_8hpp.html</filename>
    <namespace>tatami_layered</namespace>
  </compound>
  <compound kind="file">
    <name>tatami_layered.hpp</name>
    <path>tatami_layered/</path>
    <filename>tatami__layered_8hpp.html</filename>
    <namespace>tatami_layered</namespace>
  </compound>
  <compound kind="namespace">
    <name>tatami_layered</name>
    <filename>namespacetatami__layered.html</filename>
    <member kind="function">
      <type>std::shared_ptr&lt; tatami::Matrix&lt; ValueOut_, IndexOut_ &gt; &gt;</type>
      <name>convert_to_layered_sparse</name>
      <anchorfile>namespacetatami__layered.html</anchorfile>
      <anchor>aa8fecfb23f59105c5f1cadda2971afd7</anchor>
      <arglist>(const tatami::Matrix&lt; ValueIn_, IndexIn_ &gt; &amp;mat, IndexIn_ chunk_size=65536, int num_threads=1)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; tatami::Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>read_layered_sparse_from_matrix_market_text_file</name>
      <anchorfile>namespacetatami__layered.html</anchorfile>
      <anchor>a8da9e13208dd478b9fa73d3395d6e7d9</anchor>
      <arglist>(const char *filepath, Index_ chunk_size=65536, size_t buffer_size=65536)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; tatami::Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>read_layered_sparse_from_matrix_market_some_file</name>
      <anchorfile>namespacetatami__layered.html</anchorfile>
      <anchor>a48d072dc720ff2aff7ad42423af106af</anchor>
      <arglist>(const char *filepath, Index_ chunk_size=65536, size_t buffer_size=65536)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; tatami::Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>read_layered_sparse_from_matrix_market_gzip_file</name>
      <anchorfile>namespacetatami__layered.html</anchorfile>
      <anchor>aa94a50efe9d0fcefe3475eff04d2d74c</anchor>
      <arglist>(const char *filepath, Index_ chunk_size=65536, size_t buffer_size=65536)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; tatami::Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>read_layered_sparse_from_matrix_market_text_buffer</name>
      <anchorfile>namespacetatami__layered.html</anchorfile>
      <anchor>a790e143827140dcf399e7934ec8eaf7c</anchor>
      <arglist>(const unsigned char *contents, size_t length, Index_ chunk_size=65536)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; tatami::Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>read_layered_sparse_from_matrix_market_some_buffer</name>
      <anchorfile>namespacetatami__layered.html</anchorfile>
      <anchor>a01b9a0527ebb007a941c25007f6389ed</anchor>
      <arglist>(const unsigned char *contents, size_t length, Index_ chunk_size=65536, size_t buffer_size=65536)</arglist>
    </member>
    <member kind="function">
      <type>std::shared_ptr&lt; tatami::Matrix&lt; Value_, Index_ &gt; &gt;</type>
      <name>read_layered_sparse_from_matrix_market_zlib_buffer</name>
      <anchorfile>namespacetatami__layered.html</anchorfile>
      <anchor>a71bb56921e7660e8d1831e8b2d378f47</anchor>
      <arglist>(const unsigned char *contents, size_t length, Index_ chunk_size=65536, size_t buffer_size=65536)</arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>tatami helpers for creating layered matrices</title>
    <filename>index.html</filename>
    <docanchor file="index.html">md__2github_2workspace_2README</docanchor>
  </compound>
</tagfile>
