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
    filename: "main.[hash].js",
    chunkFilename: "[id].main.[hash].js",
    crossOriginLoading: "anonymous",
  },
  devtool: "inline-source-map",
  devServer: {
    compress: true,
    hot: true,
    writeToDisk: true,
  },
  plugins: [
    new webpack.HotModuleReplacementPlugin(),
    new MiniCssExtractPlugin({
      filename: "main.[hash].css",
      chunkFilename: "[id].main.[hash].css",
    }),
    new BundleTracker({
      path: __dirname,
      filename: "webpack-stats.json",
      indent: 4,
    }),
  ],
});
