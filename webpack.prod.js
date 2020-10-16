const path = require("path"),
  merge = require("webpack-merge"),
  common = require("./webpack.common.js"),
  CompressionPlugin = require("compression-webpack-plugin"),
  MiniCssExtractPlugin = require("mini-css-extract-plugin"),
  BundleTracker = require("webpack-bundle-tracker"),
  webpack = require("webpack"),
  BrotliPlugin = require("brotli-webpack-plugin");

module.exports = merge(common, {
  context: __dirname,
  mode: "production",
  output: {
    path: path.resolve(__dirname, "frontend/static/frontend/bundles/"),
    publicPath: "/static/frontend/bundles/",
    filename: "main.js",
    chunkFilename: "[id].main.js",
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
      minRatio: 0.7,
    }),
    new BrotliPlugin({
      filename: "[path].br[query]",
      algorithm: "brotliCompress",
      test: /\.(js|css|html|svg)$/,
      threshold: 10240,
      minRatio: 0.7,
    }),
  ],
});
