const path = require("path"),
  merge = require("webpack-merge"),
  common = require("./webpack.common.js"),
  CompressionPlugin = require("compression-webpack-plugin"),
  MiniCssExtractPlugin = require("mini-css-extract-plugin"),
  BundleTracker = require("webpack-bundle-tracker"),
  webpack = require("webpack");

module.exports = merge(common, {
  context: __dirname,
  mode: "production",
  output: {
    path: path.resolve(__dirname, "frontend/static/frontend/bundles/"),
    publicPath: "/static/frontend/bundles/",
    filename: "[name].[contenthash:8].js",
    chunkFilename: "[name].[contenthash:8].chunk.js",
    crossOriginLoading: "anonymous",
  },
  plugins: [
    new MiniCssExtractPlugin({
      filename: "main.css",
      chunkFilename: "[id].main.css",
    }),
    new BundleTracker({
      path: __dirname,
      filename: "webpack-stats.json",
    }),
    new webpack.optimize.SplitChunksPlugin({
      chunks: "async",
      cacheGroups: {
        vendors: {
          test: /[\\/]node_modules[\\/]/,
          priority: -10,
        },
      },
    }),
    new CompressionPlugin({
      filename: "[path].gz[query]",
      algorithm: "gzip",
      test: /\.(js|css|html)$/,
      threshold: 10240,
      minRatio: 0.8,
    }),
    new CompressionPlugin({
      filename: "[path].br[query]",
      algorithm: "brotliCompress",
      test: /\.(js|css|html|svg)$/,
      compressionOptions: {
        level: 11,
      },
      threshold: 10240,
      minRatio: 0.8,
    }),
  ],
});
