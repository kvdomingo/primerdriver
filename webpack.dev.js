const path = require("path"),
  merge = require("webpack-merge"),
  common = require("./webpack.common.js"),
  MiniCssExtractPlugin = require("mini-css-extract-plugin"),
  BundleTracker = require("webpack-bundle-tracker"),
  webpack = require("webpack");

module.exports = merge(common, {
  context: __dirname,
  mode: "development",
  output: {
    path: path.resolve(__dirname, "frontend/static/frontend/bundles/"),
    publicPath: "/static/frontend/bundles/",
    filename: "bundle.js",
    chunkFilename: "[name].chunk.js",
    crossOriginLoading: "anonymous",
  },
  devtool: "inline-source-map",
  devServer: {
    compress: true,
    hot: true,
    overlay: false,
    quiet: true,
    writeToDisk: true,
  },
  plugins: [
    new webpack.HotModuleReplacementPlugin(),
    new MiniCssExtractPlugin({
      filename: "[name].[contenthash:8].css",
      chunkFilename: "[name].[contenthash:8].chunk.css",
    }),
    new BundleTracker({
      path: __dirname,
      filename: "webpack-stats.json",
      indent: 2,
    }),
  ],
});
